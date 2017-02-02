import pyrat, glob, os, logging, copy
import numpy as np
from pyrat.viewer.Widgets import HLine, CropBoxWidget, FileselWidget, ProductContentWidget
from PyQt5 import QtCore, QtWidgets


class ESAR(pyrat.ImportWorker):
    """
    Import of DLR E-SAR SLC products
    
    This module imports E-SAR products.
    All files (i.e. SLC data and efile) must be in the same directory.
    
    :param file: The filename of the SLC file to import.
    :type array: string
    :param crop: Sets a crop area to read. Default: (0,0,0,0) = entire file
    :type array: tuple
    :author: Jens Fischer & Andreas Reigber & Michel van Kempen
    """
    gui = {'menu': 'File|Import airborne', 'entry': 'E-SAR'}
    para = [
        {'var': 'file', 'value': ''},
        {'var': 'polarisations', 'value': '*'},
        {'var': 'product', 'value': ['DCSLC', 'SLC']},
        {'var': 'crop', 'value': [0, 0, 0, 0]},
        {'var': 'bands', 'value': '*'}]
    file = ""

    def __init__(self, *args, **kwargs):
        super(ESAR, self).__init__(*args, **kwargs)
        self.name = "ESAR IMPORT"
        if len(args) == 1:
            self.file = args[0]

        if 'crop' not in self.__dict__:
            self.crop = (0, 0, 0, 0)
        else:
            if len(self.crop) != 4:
                print()
                print("    WARNING! Crop arguments must be crop=[begin_az,end_az,begin_rg,end_rg]")
                print("             ... and 0,0 can be given to read all.")
                print()

    def reader(self, *args, **kwargs):

        # if no product argument is given
        if not isinstance(self.product, str):
            if 'dcslc' in self.file:
                self.product = 'DCSLC'
                logging.info('Info: use dcslc settings')
            else:
                self.product = 'SLC'
                logging.info('Info: use slc settings')
        else:
            self.product = self.product.upper()

        if os.path.isdir(self.file):
            # polarisations aren't included in the filenames(in difference to FSAR)
            # files is a list with all files with eventually different polarisation and/or bands
            productname = 'slc' if self.product == 'SLC' else 'dc_slc'
            files = glob.glob(os.path.join(self.file, '*_ch*' + productname + '.dat'))
            files += glob.glob(os.path.join(self.file + "/*/", '*_ch*' + productname + '.dat'))

            # pick only the required files
            tmpfilelist = []
            for file in files:
                if self.bands == '*' or self.bands == crawlMetaForInfo(file, 'init.freq_band'):
                    if self.polarisations == '*' or self.polarisations == crawlMetaForInfo(file, 'init.polarization'):
                        tmpfilelist.append(file)
            files = tmpfilelist

        else:
            files = [self.file]

        # number of items to read
        number = 2 if self.product == 'SLC' else 3

        if os.path.isfile(files[0]):

            # Init array
            # ----------------------------------------------------------------
            lun = open(files[0], 'rb')
            header = np.fromfile(file=lun, dtype="int32", count=number).byteswap()
            lun.close()

            naz = header[1]
            nrg = header[0] if self.product == 'SLC' else header[2]
            block = list(self.crop)

            if block[1] == 0:
                block[1] = naz
            if block[3] == 0:
                block[3] = nrg
            daz = block[1] - block[0]
            drg = block[3] - block[2]
            offset = block[0] * nrg * np.dtype('complex64').itemsize

            array = np.empty((len(files), daz, drg), dtype='complex64')
            meta = {}
            meta['CH_pol'] = [' '] * len(files)
        # ----------------------------------------------------------------

        else:
            logging.error("File not found: " + ESAR.file)
            return

        for k, file in enumerate(files):
            logging.info("Found " + file)

            # READ DATA
            if os.path.isfile(file):

                # Read array
                # ---------------------------------------------------------------------------------------
                lun = open(file, 'rb')
                header = np.fromfile(file=lun, dtype="int32", count=number).byteswap()
                lun.seek(offset, 1)

                npix = daz * nrg
                foo = np.fromfile(file=lun, dtype="complex64", count=npix).byteswap().reshape(daz, nrg)
                array[k, ...] = foo[:, block[2]:block[3]]
                lun.close()
                # ---------------------------------------------------------------------------------------

            else:
                logging.error("File not found: " + file)

            # READ EFILE - METADATA

            fl = os.path.split(file)
            efile = fl[0] + '/e' + fl[1][1:].replace('_slc.dat' if self.product == 'SLC' else '_dcslc.dat', '.txt')
            if os.path.isfile(efile):

                lun = open(efile)
                etext = lun.readlines()
                lun.close()

                meta['sensor'] = 'DLR E-SAR'
                meta['CH_pol'][k] = next(line for line in etext if 'init.polarization' in line)[-3:-1]

                if k == 0:
                    meta['c0'] = float(next(line for line in etext if 'init.speed_of_light' in line)[-15:])
                    meta['rd'] = float(next(line for line in etext if 'init.range_delay' in line)[-15:]) * 1e-6
                    meta['rsf'] = float(next(line for line in etext if 'init.range_sampling_rate' in line)[-15:]) * 1e6
                    meta['nrg'] = drg
                    meta['naz'] = daz
                    meta['lam'] = float(next(line for line in etext if 'init.wavelength' in line)[-15:])
                    meta['band'] = next(line for line in etext if 'init.freq_band' in line)[-15:].strip()
                    if meta['band'] == 'X':
                        meta['antdir'] = +1
                    else:
                        meta['antdir'] = -1
                    meta['v0'] = float(next(line for line in etext if 'init.forw_velocity' in line)[-15:])
                    meta['bw'] = float(next(line for line in etext if 'init.chirp_bandwidth' in line)[-15:]) * 1e6

                    meta['rd'] += block[2] / meta['rsf']
                    az_presum_post_slc = float(next(line for line in etext if 'init.az_presum_post_slc' in line)[-15:])
                    az_presum_pre = float(next(line for line in etext if 'init.az_presum_pre' in line)[-15:])
                    meta['prf'] = float(
                        next(line for line in etext if 'init.prf' in line)[-15:]) / az_presum_post_slc / az_presum_pre

                    res_az_ml = float(next(line for line in etext if 'init.reso_azimuth' in line)[-15:])
                    nlooks = float(next(line for line in etext if 'init.looks' in line)[-15:])

                    bw_rg = float(next(line for line in etext if 'init.chirp_bandwidth' in line)[-15:]) * 1e6

                    meta['pre_az'] = az_presum_pre * az_presum_post_slc
                    meta['res_az'] = res_az_ml / ((nlooks + 1) / 2.)
                    meta['res_rg'] = 1.33 * (meta['c0'] / 2.) / bw_rg
            else:
                logging.error("Metadata file not found: " + efile)

        return array, meta

    @classmethod
    def guirun(cls, viewer):
        para_backup = copy.deepcopy(cls.para)  # keep a deep copy of the default parameters
        wid = EsarImportWidget(dir=cls.para[[par['var'] for par in cls.para].index('file')]['value'])
        wid.update()
        res = wid.exec_()
        if res == 1:
            plugin = cls(filename=wid.dir, product=wid.product, bands=wid.band, polarisations=wid.polar, crop=wid.crop)
            viewer.statusBar.setMessage(message=plugin.name + ' running', colour='R')
            plugin.run()
            del plugin
            viewer.statusBar.setMessage(message='Ready', colour='G')
            viewer.updateViewer()


class EsarImportWidget(QtWidgets.QDialog):
    def __init__(self, parent=None, dir=None):
        super(EsarImportWidget, self).__init__(parent)
        self.setWindowTitle("ESAR import")
        mainlayout = QtWidgets.QVBoxLayout(self)

        self.dirwidget = FileselWidget(title='ESAR product directory (RGI-SR)', type='opendir')
        self.dirwidget.setvalue(dir)
        mainlayout.addWidget(self.dirwidget)
        mainlayout.addWidget(HLine())
        self.productwidget = ProductContentWidget(products=["SLC", "DCSLC"])
        mainlayout.addWidget(self.productwidget)
        mainlayout.addWidget(HLine())
        self.cropwidget = CropBoxWidget(title='Select crop (0=maximum)')
        mainlayout.addWidget(self.cropwidget)

        self.buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel,
                                              QtCore.Qt.Horizontal, self)
        mainlayout.addWidget(self.buttons)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

        self.dirwidget.text.textChanged.connect(lambda: self.update(mode=0))  # update all

        self.bandupdate = lambda: self.update(mode=2)  # update band
        self.productupdate = lambda: self.update(mode=1)
        self.productwidget.product.currentIndexChanged.connect(self.productupdate)  # update product
        self.productwidget.band.currentIndexChanged.connect(self.bandupdate)

    def update(self, mode=0):
        self.dir = str(self.dirwidget.getvalue())
        ESAR.file = self.dir
        self.product = self.productwidget.getvalue(0)
        self.polar = self.productwidget.getvalue(2)

        if os.path.isdir(self.dir):  # read from several files
            allfiles = glob.glob(os.path.join(self.dir, '*slc*' + '*.dat'))
        elif os.path.isfile(self.dir):  # read from a single file
            allfiles = [self.dir]
        else:  # nothing selected
            allfiles = []

        # allpols = list(set([os.path.basename(slc).split('_')[1][2:].upper() for slc in allfiles]))
        allpols = list(set([crawlMetaForInfo(slc, 'init.polarization') for slc in allfiles]))
        allbands = list(set([crawlMetaForInfo(slc, 'init.freq_band') for slc in allfiles]))
        allproducts = list(set(['SLC' if 'slc' in slc else 'DCSLC' for slc in allfiles]))

        # the range settings for the crop widget are missing (the max value is 99999...)
        # take a look on F-SAR; where are the dimensions of an image (.dat) from E-SAR?
        # self.cropwidget.setrange([[0, naz], [0, naz], [0, nrg], [0, nrg]])
        # self.cropwidget.setvalues([0, naz, 0, nrg])

        if mode == 0 or mode == 1:
            self.productwidget.band.currentIndexChanged.disconnect(self.bandupdate)
            self.productwidget.product.currentIndexChanged.disconnect(self.productupdate)
            self.productwidget.updatebands(allbands)
            self.productwidget.updatepolar(allpols)
            self.productwidget.updateproducts(allproducts)
            self.productwidget.setvalue(2, self.polar)
            self.productwidget.band.currentIndexChanged.connect(self.bandupdate)
            self.productwidget.product.currentIndexChanged.connect(self.productupdate)

        elif mode == 2:
            self.productwidget.updatepolar(allpols)
            self.productwidget.setvalue(2, self.polar)

    def accept(self):
        self.product = self.productwidget.getvalue(0)
        self.polar = self.productwidget.getvalue(2)
        self.crop = self.cropwidget.getvalues()
        self.band = self.productwidget.getvalue(1)
        super(EsarImportWidget, self).accept()


def crawlMetaForInfo(path, nametag):
    """
    reads the value with the nametag from the meta file, which belongs to the *.dat file and returns it as string
    only written for ESAR-metadata (.txt) files,
    the meta file should be in the same dir...
    :param path: absolute path from the .dat file
    :param nametag: name from the attribute (for example 'init.polarization')
    """
    try:
        # */i*_slc.dat --> */e*.txt
        fl = os.path.split(path)
        efile = fl[0] + '/e' + fl[1][1:].replace('_slc.dat', '.txt')

        # are the data files in a dir which names 'RGI-SR'? - solution for a specific case
        if os.path.split(efile)[0][-6:] == 'RGI-SR':
            # then change the directory to 'RGI-RDP'
            efile = efile.replace('RGI-SR', 'RGI-RDP')

        # get meta data
        lun = open(efile, 'r')
        etext = lun.readlines()
        lun.close()

        # find the line in meta
        data = next(line for line in etext if nametag in line)
        if data == None:
            logging.error("Warning: no " + nametag + " found in metadata file")
            return

        data = data[data.find(nametag) + len(nametag):]
        data = data.strip()

        return data

    except FileNotFoundError:
        # logging.error("Warning: perhaps the format isn't completely supported, some features wont function ")
        return None


@pyrat.docstringfrom(ESAR)
def esar(*args, **kwargs):
    return ESAR(*args, **kwargs).run(*args, **kwargs)


class ESAR_track(pyrat.ImportWorker):
    """
    Import of DLR E-SAR REF/REALTRACK products

    :author: Andreas Reigber
    :param file: The filename of the DCSLC file to import.
    :type array: string
    :param setlen: Interpolate to the given azimuth length.
    :type array: int
    :param crop: Sets a crop area to read. Default: (0,0,0,0) = entire file
    :type crop: tuple

    """

    def __init__(self, *args, **kwargs):
        super(ESAR_track, self).__init__(*args, **kwargs)
        self.name = "ESAR REAL/REF TRACK IMPORT"
        self.block = 'T'
        if len(args) == 1:
            self.file = args[0]
        if 'crop' not in self.__dict__:  self.crop = (0, 0, 0, 0)

    def reader(self, *args, **kwargs):
        meta = {}
        meta['sensor'] = 'DLR E-SAR'
        array = False

        if os.path.isfile(self.file):
            lun = open(self.file, 'rb')
            nry = np.fromfile(file=lun, dtype="int32", count=1).byteswap()[0]
            time = np.fromfile(file=lun, dtype="float64", count=nry).byteswap()
            trk1 = np.fromfile(file=lun, dtype="float64", count=nry * 3).byteswap().reshape(3, nry)
            trk2 = np.fromfile(file=lun, dtype="float64", count=nry * 3).byteswap().reshape(3, nry)
            trk3 = np.fromfile(file=lun, dtype="float64", count=nry * 3).byteswap().reshape(3, nry)
            lun.close()

            if hasattr(self, 'setlen'):
                naz = self.setlen
                track = np.empty((naz, 7), dtype='float64')
                time = np.interp(np.arange(naz, dtype='float64') / naz * nry, np.arange(nry, dtype='float64'), time)
                for l in range(3):
                    track[:, l + 0] = np.interp(np.arange(naz, dtype='float64') / naz * nry,
                                                np.arange(nry, dtype='float64'), trk2[l, :])
                    track[:, l + 3] = np.interp(np.arange(naz, dtype='float64') / naz * nry,
                                                np.arange(nry, dtype='float64'), trk1[l, :])
                track[:, 6] = time
            else:
                naz = nry
                track = np.empty((naz, 7), dtype='float64')
                for l in range(3):
                    track[:, l + 0] = trk2[l, :]
                    track[:, l + 3] = trk1[l, :]
                track[:, 6] = time

            block = list(self.crop)
            if block[1] == 0: block[1] = naz
            return track[block[0]:block[1], :], meta
        else:
            logging.error("Track file not found (setting track to None): " + self.file)
            return None, None


@pyrat.docstringfrom(ESAR_track)
def esar_track(*args, **kwargs):
    return ESAR_track(*args, **kwargs).run(*args, **kwargs)
