import pyrat
import glob, os
import logging
import copy
import numpy as np
from PyQt5 import QtCore, QtWidgets

# from pyrat.load import RatFile
from pyrat.lib.ste import RatFile, Xml2Py
from pyrat.viewer.Dialogs import FlexFilesel
from pyrat.viewer.Widgets import HLine, CropBoxWidget, BoolWidget, FileselWidget, ProductContentWidget


class FSAR(pyrat.ImportWorker):
    """
    Import of DLR F-SAR SLC product. This class loads the SLC data of one or several bands and / or one
    or several polarisations, together with their meta data into a new PyRAT layer(s).

    :param dir: The F-SAR product directory.
    :type dir: str
    :param bands: Load only this band. '*' to load all bands. Default='*'
    :type band: string
    :param polarisations: Load only this polarisation. '*' to load all bands. Default='*'
    :type polarisation: string
    :param product: Selects the product component to import. Default='RGI-SLC'
    :type product: string
    :param suffix: An optional suffix appended to the INF directory (i.e. INF_<suffix>)
    :type suffix: string
    :param crop: A crop region / subset to import (az_start, az_end, rg_start, rg_end)
    :type crop: tuple
    :param sym: PolSAR symmetrisation. If set, HV and VH are averaged on import. Default=False
    :type sym: bool
    :author: Andreas Reigber
    """
    gui = {'menu': 'File|Import airborne', 'entry': 'F-SAR'}
    para = [
        {'var': 'dir', 'value': ''},
        {'var': 'bands', 'value': '*'},
        {'var': 'polarisations', 'value': '*'},
        {'var': 'product', 'value': 'RGI-SLC'},
        {'var': 'suffix', 'value': ''},
        {'var': 'crop', 'value': [0, 0, 0, 0]},
        {'var': 'sym', 'value': False, 'type': 'bool', 'text': 'Cross-polar symmetrisation'},
        {'var': 'mask', 'type': bool, 'value': False}]

    def __init__(self, *args, **kwargs):
        super(FSAR, self).__init__(*args, **kwargs)
        self.name = "FSAR SLC IMPORT"
        if len(args) == 1:
            self.dir = args[0]

    def reader(self, *args, **kwargs):

        if self.product == 'RGI-SLC':
            head = 'slc'
            src = ('RGI', 'RGI-SR')
        if self.product == 'RGI-AMP':
            head = 'amp'
            src = ('RGI', 'RGI-SR')
        if self.product == 'INF-SLC':
            head = 'slc_coreg_offset'
            src = ('INF' + ('_' + self.suffix if len(self.suffix) > 0 else ''), 'INF-SR')
        if self.product == 'INF-CIR':
            head = 'slcpol'
            src = ('INF', 'INF-SR')

        if self.polarisations == '*':
            self.polarisations = '??'
        if self.product == 'INF-CIR':
            files = glob.glob(os.path.join(self.dir, src[0], src[1],
                                           head + '*' + self.bands.upper() + self.polarisations.lower() + '*_c' + str(
                                               self.track) + '_s' + str(self.subaperture) + '_*coreg.rat'))
        else:
            files = glob.glob(os.path.join(self.dir, src[0], src[1],
                                           head + '*' + self.bands.upper() + self.polarisations.lower() + '_*.rat'))

            if self.product == 'INF-SLC':
                if len(files) == 0:
                    head = 'slc_coreg'
                    files = glob.glob(os.path.join(self.dir, src[0], src[1], head + '*' + self.bands.upper()
                                                   + self.polarisations.lower() + '_*.rat'))

            bands = list(set([os.path.basename(slc).split('_')[-2][0] for slc in files]))
            pols = list(set([os.path.basename(slc).split('_')[-2][1:3] for slc in files]))

        if self.product != 'INF-SLC':
            bands = list(set([os.path.basename(slc).split('_')[2][0] for slc in files]))
            pols = list(set([os.path.basename(slc).split('_')[2][1:3] for slc in files]))


        array = []
        meta = []
        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band + '-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_' + band in f]
                fil = RatFile(bandfiles[0])
                if self.mask is True:
                    maskfile = glob.glob(os.path.join(self.dir, src[0], src[1], 'mask*' + band.upper() + '*.rat'))
                    msk = RatFile(maskfile[0])
                naz = fil.shape[0]
                nrg = fil.shape[1]
                block = list(self.crop)
                if block[1] == 0 or block[1] > naz:
                    block[1] = naz
                if block[3] == 0 or block[3] > nrg:
                    block[3] = nrg
                daz = block[1] - block[0]
                drg = block[3] - block[2]

                barr = np.empty((len(bandfiles), daz, drg), dtype='complex64')
                bmeta = {}
                bmeta['sensor'] = 'DLR F-SAR'
                bmeta['band'] = band
                bmeta['CH_pol'] = [' '] * len(bandfiles)
                for k, f in enumerate(bandfiles):
                    logging.info("Found " + f)
                    barr[k, ...] = RatFile(f).read(block=block)
                    if self.mask is True:
                        mask = msk.read(block=block)
                        # barr[k, ...] *= mask
                    if self.product == 'RGI-SLC':
                        ppfile = f.replace('RGI-SR', 'RGI-RDP').replace('slc_', 'pp_').replace('.rat', '.xml')
                    if self.product == 'RGI-AMP':
                        ppfile = f.replace('RGI-SR', 'RGI-RDP').replace('amp_', 'pp_').replace('.rat', '.xml')
                    if self.product == 'INF-SLC':
                        ppname = 'pp_' + '_'.join(os.path.basename(f).split('_')[-3:]).replace('.rat', '.xml')
                        ppfile = os.path.join(self.dir, src[0], 'INF-RDP', ppname)
                    if self.product == 'INF-CIR':
                        ppname = 'ppgeo_csar_' + '_'.join(os.path.basename(f).split('_')[1:4]) + '.xml'
                        ppfile = os.path.join(self.dir, 'GTC', 'GTC-RDP', ppname)
                    pp = Xml2Py(ppfile)
                    bmeta['CH_pol'][k] = pp.polarisation

                if self.sym is True and barr.ndim == 3 and barr.shape[0] == 4:
                    pol = bmeta['CH_pol']
                    idx_hh = pol.index('HH')
                    idx_vv = pol.index('VV')
                    idx_hv = pol.index('HV')
                    idx_vh = pol.index('VH')
                    barr[idx_hv, ...] =  (barr[idx_hv, ...] + barr[idx_vh, ...]) / np.sqrt(2)
                    barr = np.delete(barr, idx_vh, axis=0)
                    bmeta['CH_pol'][idx_hv] = 'XX'
                    bmeta['CH_pol'].remove('VH')

                bmeta['prf'] = pp.prf
                bmeta['c0'] = pp.c0
                bmeta['rd'] = pp.rd
                bmeta['rsf'] = pp.rsf
                bmeta['nrg'] = drg
                bmeta['naz'] = daz
                bmeta['lam'] = pp['lambda']
                bmeta['band'] = pp.band
                bmeta['antdir'] = pp.antdir
                bmeta['v0'] = pp.v0
                bmeta['bw'] = pp.cbw
                bmeta['ps_rg'] = pp.ps_rg
                bmeta['ps_az'] = pp.ps_az
                bmeta['rd'] += block[2] / bmeta['rsf']
                bmeta['h0'] = pp.h0
                bmeta['pre_az'] = pp.pre_az
                bmeta['terrain'] = pp.terrain

                if self.mask is True:
                    array.append([barr, mask])
                else:
                    array.append(barr)
                meta.append(bmeta)

        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], meta[0]
        else:
            return array, meta

    @classmethod
    def guirun(cls, viewer):
        para_backup = copy.deepcopy(cls.para)  # keep a deep copy of the default parameters
        wid = FsarImportWidget()
        wid.update()
        res = wid.exec_()
        if res == 1:
            plugin = cls(dir=wid.dir, product=wid.product, bands=wid.bands, polarisations=wid.polar, crop=wid.crop, sym=wid.sym)
            viewer.statusBar.setMessage(message=plugin.name + ' running', colour='R')
            plugin.run()
            del plugin
            viewer.statusBar.setMessage(message='Ready', colour='G')
            viewer.updateViewer()


class FsarImportWidget(QtWidgets.QDialog):
    def __init__(self, parent=None, dir=None):
        super(FsarImportWidget, self).__init__(parent)
        self.setWindowTitle("FSAR import")
        mainlayout = QtWidgets.QVBoxLayout(self)

        self.dirwidget = FileselWidget(title='FSAR product dir', type='opendir')
        self.dirwidget.setvalue(dir)
        mainlayout.addWidget(self.dirwidget)
        mainlayout.addWidget(HLine())
        self.productwidget = ProductContentWidget(products=["RGI-SLC", "RGI-AMP", "INF-SLC"])
        mainlayout.addWidget(self.productwidget)
        mainlayout.addWidget(HLine())
        self.cropwidget = CropBoxWidget(title='Select crop (0=maximum)')
        mainlayout.addWidget(self.cropwidget)
        mainlayout.addWidget(HLine())
        self.symwidget = BoolWidget(text="Cross-polar symmetrisation")
        mainlayout.addWidget(self.symwidget)

        self.buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel,
                                                  QtCore.Qt.Horizontal, self)
        mainlayout.addWidget(self.buttons)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

        self.dirwidget.text.textChanged.connect(lambda: self.update(mode=0))  # update all

        self.bandupdate = lambda: self.update(mode=2)  # update band
        self.productwidget.product.currentIndexChanged.connect(lambda: self.update(mode=1))  # update product
        self.productwidget.band.currentIndexChanged.connect(self.bandupdate)

    def update(self, mode=0):
        self.dir = str(self.dirwidget.getvalue())
        self.product = self.productwidget.getvalue(0)
        self.bands = self.productwidget.getvalue(1)
        self.polar = self.productwidget.getvalue(2)

        if self.product == 'RGI-SLC':
            head = 'slc'
            src = ('RGI', 'RGI-SR')
            code_pos = 2
        if self.product == 'RGI-AMP':
            head = 'amp'
            src = ('RGI', 'RGI-SR')
            code_pos = 2
        if self.product == 'INF-SLC':
            head = 'slc_coreg'
            src = ('INF', 'INF-SR')
            code_pos = 4

        files = glob.glob(
            os.path.join(self.dir, src[0], src[1], head + '*' + self.bands.upper() + self.polar.lower() + '*.rat'))

        if mode == 0:
            allfiles = glob.glob(os.path.join(self.dir, src[0], src[1], head + '*.rat'))
            self.bands = '*'
            self.polar = '*'
        if mode == 1:
            allfiles = glob.glob(os.path.join(self.dir, src[0], src[1], head + '*.rat'))
        if mode == 2:
            allfiles = glob.glob(os.path.join(self.dir, src[0], src[1], head + '*' + self.bands.upper() + '*.rat'))

        # allfiles = glob.glob(os.path.join(self.dir,'RGI','RGI-SR',head+'*.rat'))
        allbands = list(set([os.path.basename(slc).split('_')[code_pos][0] for slc in allfiles]))
        allpols = list(set([os.path.basename(slc).split('_')[code_pos][1:3].upper() for slc in allfiles]))

        nrg = 0
        naz = 0
        for filename in files:
            lun = RatFile(filename)
            nrg = max(nrg, lun.shape[1])
            naz = max(naz, lun.shape[0])
        self.cropwidget.setrange([[0, naz], [0, naz], [0, nrg], [0, nrg]])
        self.cropwidget.setvalues([0, naz, 0, nrg])

        if mode == 0 or mode == 1:
            self.productwidget.band.currentIndexChanged.disconnect(self.bandupdate)
            self.productwidget.updatepolar(allpols)
            self.productwidget.updatebands(allbands)
            self.productwidget.setvalue(1, self.bands)
            self.productwidget.setvalue(2, self.polar)
            self.productwidget.band.currentIndexChanged.connect(self.bandupdate)

        elif mode == 2:
            self.productwidget.updatepolar(allpols)
            self.productwidget.setvalue(2, self.polar)

    def accept(self):
        self.product = self.productwidget.getvalue(0)
        self.bands = self.productwidget.getvalue(1)
        self.polar = self.productwidget.getvalue(2)
        self.crop = self.cropwidget.getvalues()
        self.sym = self.symwidget.getvalue()
        super(FsarImportWidget, self).accept()


@pyrat.docstringfrom(FSAR)
def fsar(*args, **kwargs):
    return FSAR(*args, **kwargs).run(*args, **kwargs)


class FSAR_dem(pyrat.ImportWorker):
    """
    Import of DLR F-SAR SLC product

    :param dir: The F-SAR product directory.
    :type dir: str
    :param match: A matching string to select subset of files
    :type match: string
    :param crop: A crop region / subset to import (az_start, az_end, rg_start, rg_end)
    :type crop: tuple
    :author: Andreas Reigber
    """

    para = [
        {'var': 'dir', 'value': ''},
        {'var': 'bands', 'value': '*'},
        {'var': 'polarisations', 'value': '*'},
        {'var': 'product', 'value': 'RGI-SLC'},
        {'var': 'crop', 'value': [0, 0, 0, 0]}
    ]

    def __init__(self, *args, **kwargs):
        super(FSAR_dem, self).__init__(*args, **kwargs)
        self.name = "FSAR DEM IMPORT"
        if len(args) == 1:
            self.dir = args[0]

    def reader(self, *args, **kwargs):

        head = 'slantdem_full'
        src = ('RGI', 'RGI-AUX')

        files = glob.glob(os.path.join(self.dir, src[0], src[1], head + '*' + self.bands.upper() + '*.rat'))
        bands = list(set([os.path.basename(slc).split('_')[2][0] for slc in files]))
        pols = list(set([os.path.basename(slc).split('_')[2][1:3] for slc in files]))

        array = []
        meta = []
        meta.append({})
        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band + '-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_' + band in f]
                fil = RatFile(bandfiles[0])
                naz = fil.shape[0]
                nrg = fil.shape[1]
                block = list(self.crop)
                if block[1] == 0:
                    block[1] = naz
                if block[3] == 0:
                    block[3] = nrg
                daz = block[1] - block[0]
                drg = block[3] - block[2]

                barr = np.empty((len(bandfiles), daz, drg), dtype='float32')
                for k, f in enumerate(bandfiles):
                    logging.info("Found " + f)
                    barr[k, ...] = RatFile(f).read(block=block)
                array.append(barr)

        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], meta[0]
        else:
            return array, meta


@pyrat.docstringfrom(FSAR_dem)
def fsar_dem(*args, **kwargs):
    return FSAR_dem(*args, **kwargs).run(*args, **kwargs)


class FSAR_offnadir(pyrat.ImportWorker):
    """
    Import of DLR F-SAR SLC product

    :param dir: The F-SAR product directory.
    :type dir: str
    :param match: A matching string to select subset of files
    :type match: string
    :param crop: A crop region / subset to import (az_start, az_end, rg_start, rg_end)
    :type crop: tuple
    :author: Andreas Reigber
    """

    para = [
        {'var': 'dir', 'value': ''},
        {'var': 'bands', 'value': '*'},
        {'var': 'polarisations', 'value': '*'},
        {'var': 'product', 'value': 'RGI-SLC'},
        {'var': 'crop', 'value': [0, 0, 0, 0]}
    ]

    def __init__(self, *args, **kwargs):
        super(FSAR_offnadir, self).__init__(*args, **kwargs)
        self.name = "FSAR OFFNADIR IMPORT"
        if len(args) == 1:
            self.dir = args[0]

    def reader(self, *args, **kwargs):

        head = 'offnadir'
        src = ('RGI', 'RGI-AUX')

        files = glob.glob(os.path.join(self.dir, src[0], src[1], head + '*' + self.bands.upper() + '*.rat'))
        bands = list(set([os.path.basename(slc).split('_')[2][0] for slc in files]))
        pols = list(set([os.path.basename(slc).split('_')[2][1:3] for slc in files]))

        array = []
        meta = []
        meta.append({})
        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band + '-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_' + band in f]
                fil = RatFile(bandfiles[0])
                naz = fil.shape[0]
                nrg = fil.shape[1]
                block = list(self.crop)
                if block[1] == 0:
                    block[1] = naz
                if block[3] == 0:
                    block[3] = nrg
                daz = block[1] - block[0]
                drg = block[3] - block[2]

                barr = np.empty((len(bandfiles), daz, drg), dtype='float32')
                for k, f in enumerate(bandfiles):
                    logging.info("Found " + f)
                    barr[k, ...] = RatFile(f).read(block=block)
                array.append(np.squeeze(barr))
        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], meta[0]
        else:
            return array, meta


@pyrat.docstringfrom(FSAR_offnadir)
def fsar_offnadir(*args, **kwargs):
    return FSAR_offnadir(*args, **kwargs).run(*args, **kwargs)


class FSAR_phadem(pyrat.ImportWorker):
    """
    Import of DLR F-SAR DEM PHASE product

    :param dir: The F-SAR product directory.
    :type dir: str
    :param match: A matching string to select subset of files
    :type match: string
    :param suffix: An optional suffix appended to the INF directory (i.e. INF_<suffix>)
    :type suffix: string
    :param crop: A crop region / subset to import (az_start, az_end, rg_start, rg_end)
    :type crop: tuple
    :author: Andreas Reigber
    """

    para = [
        {'var': 'dir', 'value': ''},
        {'var': 'bands', 'value': '*'},
        {'var': 'product', 'value': 'RGI-SLC'},
        {'var': 'suffix', 'value': ''},
        {'var': 'crop', 'value': [0, 0, 0, 0]}
    ]

    def __init__(self, *args, **kwargs):
        super(FSAR_phadem, self).__init__(*args, **kwargs)
        self.name = "FSAR DEM PHASE IMPORT"
        if len(args) == 1:
            self.dir = args[0]

    def reader(self, *args, **kwargs):

        head = 'pha_dem'
        src = ('INF' + ('_' + self.suffix if len(self.suffix) > 0 else ''), 'INF-SR')

        files = glob.glob(os.path.join(self.dir, src[0], src[1], head + '*' + self.bands.upper() + '*.rat'))
        bands = list(set([os.path.basename(slc).split('_')[2][0] for slc in files]))

        array = []
        meta = [{}]
        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band + '-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_' + band in f]
                fil = RatFile(bandfiles[0])
                naz = fil.shape[0]
                nrg = fil.shape[1]
                block = list(self.crop)
                if block[1] == 0:
                    block[1] = naz
                if block[3] == 0:
                    block[3] = nrg
                daz = block[1] - block[0]
                drg = block[3] - block[2]

                barr = np.empty((len(bandfiles), daz, drg), dtype='float32')
                for k, f in enumerate(bandfiles):
                    logging.info("Found " + f)
                    barr[k, ...] = RatFile(f).read(block=block)
                array.append(barr)

        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], meta[0]
        else:
            return array, meta


@pyrat.docstringfrom(FSAR_phadem)
def fsar_phadem(*args, **kwargs):
    return FSAR_phadem(*args, **kwargs).run(*args, **kwargs)


class FSAR_kz(pyrat.ImportWorker):
    """
    Import of DLR F-SAR KZ interferometric product

    :param dir: The F-SAR product directory.
    :type dir: str
    :param match: A matching string to select subset of files
    :type match: string
    :param suffix: An optional suffix appended to the INF directory (i.e. INF_<suffix>)
    :type suffix: string
    :param crop: A crop region / subset to import (az_start, az_end, rg_start, rg_end)
    :type crop: tuple
    :author: Andreas Reigber
    """

    para = [
        {'var': 'dir', 'value': ''},
        {'var': 'bands', 'value': '*'},
        {'var': 'polarisations', 'value': '*'},
        {'var': 'suffix', 'value': ''},
        {'var': 'sym', 'value': False, 'type': 'bool', 'text': 'Cross-polar symmetrisation'},
        {'var': 'crop', 'value': [0, 0, 0, 0]}
    ]

    def __init__(self, *args, **kwargs):
        super(FSAR_kz, self).__init__(*args, **kwargs)
        self.name = "FSAR KZ IMPORT"
        if len(args) == 1:
            self.dir = args[0]

    def reader(self, *args, **kwargs):

        head = 'kz'
        src = ('INF' + ('_' + self.suffix if len(self.suffix) > 0 else ''), 'INF-SR')

        files = glob.glob(os.path.join(self.dir, src[0], src[1], head + '*'
                                       + self.bands.upper() + self.polarisations.lower() + '*.rat'))
        bands = list(set([os.path.basename(slc).split('_')[3][0] for slc in files]))
        pols = list(set([os.path.basename(slc).split('_')[3][1:3] for slc in files]))

        array = []
        meta = []
        meta.append({})

        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band + '-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_' + band in f]
                if self.sym is True and len(bandfiles) == 4:           # Remove HV channel (data was symmetrised)
                    for slc in bandfiles:
                        if os.path.basename(slc).split('_')[3][1:3] == 'vh':
                            bandfiles.remove(slc)

                fil = RatFile(bandfiles[0])
                naz = fil.shape[0]
                nrg = fil.shape[1]
                block = list(self.crop)
                if block[1] == 0:
                    block[1] = naz
                if block[3] == 0:
                    block[3] = nrg
                daz = block[1] - block[0]
                drg = block[3] - block[2]

                barr = np.empty((len(bandfiles), daz, drg), dtype='float32')
                for k, f in enumerate(bandfiles):
                    logging.info("Found " + f)
                    barr[k, ...] = RatFile(f).read(block=block)
                array.append(barr)

        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], meta[0]
        else:
            return array, meta


@pyrat.docstringfrom(FSAR_dem)
def fsar_kz(*args, **kwargs):
    return FSAR_kz(*args, **kwargs).run(*args, **kwargs)


class FSAR_track(pyrat.ImportWorker):
    """
    Import of DLR E-SAR REF/REALTRACK products

    :param dir: The F-SAR product directory.
    :type dir: str
    :param match: A matching string to select subset of files
    :type match: string
    :param crop: A crop region / subset to import (az_start, az_end, rg_start, rg_end)
    :type crop: tuple
    :author: Andreas Reigber
    """

    para = [
        {'var': 'dir', 'value': ''},
        {'var': 'bands', 'value': '*'},
        {'var': 'polarisations', 'value': '*'},
        {'var': 'product', 'value': 'RGI-SLC'},
        {'var': 'crop', 'value': [0, 0, 0, 0]}]

    def __init__(self, *args, **kwargs):
        super(FSAR_track, self).__init__(*args, **kwargs)
        self.name = "FSAR REAL/REF TRACK IMPORT"
        self.block = 'T'
        if len(args) == 1:
            self.dir = args[0]

    def reader(self, *args, **kwargs):
        sys = pyrat.data.getAnnotation()
        pre_az = sys['pre_az']  # pressuming in azimuth factor
        if self.product == 'RGI-SLC':
            # head = 'reftr_sar_resa'
            head = '_sar_resa'
            src = ('RGI', 'RGI-TRACK')
        if self.product == 'RGI-AMP':
            # head = 'reftr_sar_resa'
            head = '_sar_resa'
            src = ('RGI', 'RGI-TRACK')
        if self.product == 'INF-SLC':
            # head = 'reftr_sar_resa'
            head = '_sar_resa'
            src = ('INF', 'INF-TRACK')
        if self.product == 'INF-CIR':
            # head = 'reftr_sar_resa'
            head = 'track_loc'
            src = ('GTC', 'GTC-AUX')

        files = glob.glob(os.path.join(self.dir, src[0], src[1],
                                       '*' + head + '*' + self.bands.upper() + self.polarisations.lower() + '*.rat'))
        if self.product == 'INF-CIR':
            bands = list(set([os.path.basename(slc).split('_')[3][0] for slc in files]))
            pols = list(set([os.path.basename(slc).split('_')[3][1:3] for slc in files]))
            tracks = list(set([os.path.basename(slc).split('_')[0][0:6] for slc in
                               files]))  # contains the path to the 2 track files (reference and real)
        else:
            bands = list(set([os.path.basename(slc).split('_')[4][0] for slc in files]))
            pols = list(set([os.path.basename(slc).split('_')[4][1:3] for slc in files]))
            tracks = list(set([os.path.basename(slc).split('_')[0][0:6] for slc in
                               files]))  # contains the path to the 2 track files (reference and real)

        array = []
        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band + '-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_' + band in f]  # list of files from the same band
                fil = RatFile(bandfiles[0])
                naz = fil.shape[-2]
                block = list(self.crop)
                block[2] = 0
                block[3] = 7
                if block[1] == 0:
                    block[1] = naz
                daz = block[1] - block[0]
                drg = block[3] - block[2]

                barr = np.empty((len(bandfiles) / 2, sys['pre_az'] * daz, 7))
                for k, pol in enumerate(pols):
                    polfiles = [f for f in bandfiles if pol + '_' in f]
                    for i, f in enumerate(polfiles):
                        logging.info("Found " + f)
                        if 'reftr' in f:
                            barr[k, :, 0:3] = RatFile(f).read(block=(pre_az * block[0], 1, pre_az * daz, 3)).T
                            barr[k, :, 6] = RatFile(f).read(block=(pre_az * block[0], 0, pre_az * daz, 1))
                        elif 'track' in f:
                            barr[k, :, 3:6] = RatFile(f).read(block=(pre_az * block[0], 1, pre_az * daz, 3)).T
                            if self.product == 'INF-SLC':  # read multisquit if existing
                                dir = f.split('INF/INF-TRACK/track_sar_resa')[0] + 'INF/INF-AUX/'
                                master = f.split('_')[-3]
                                band = f.split('_')[-2][0]
                                ms_files = glob.glob(dir + 'baseline_error_*' + master + '*' + band + '*.rat')
                                for msf in ms_files:
                                    if 'constLin' in msf:
                                        logging.info('Mutisquint const/linear update found!')
                                        ms_corr = RatFile(msf).read()
                                        x = barr[k, :, 3] - barr[k, 0, 3]
                                        barr[k, :, 4] -= (ms_corr[0] + x * ms_corr[1])
                                        barr[k, :, 5] -= (ms_corr[2] + x * ms_corr[3])
                                    else:
                                        logging.info('Mutisquint baseline correction found!')
                                        ms_corr = RatFile(msf).read()[..., block[0]:block[0] + daz]
                                        dy = np.sum(ms_corr[1, ...], axis=0)
                                        dz = np.sum(ms_corr[2, ...], axis=0)
                                        barr[k, :, 4] += np.resize(dy, pre_az * daz)
                                        barr[k, :, 5] += np.resize(dz, pre_az * daz)

                array.append(barr[:, 0:pre_az * daz:pre_az, :])

        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], None
        else:
            return array, None


@pyrat.docstringfrom(FSAR_track)
def fsar_track(*args, **kwargs):
    return FSAR_track(*args, **kwargs).run(*args, **kwargs)
