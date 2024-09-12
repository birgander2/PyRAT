"""
Import UAVSAR SLC's

Author: Gustavo Daniel Martín del Campo Becerra
E-Mail: Gustavo.MartindelCampoBecerra@dlr.de
"""

import pyrat
import glob, os
import logging
import copy
import numpy as np
from PyQt5 import QtCore, QtWidgets
from pyrat.tools import bcolors
from pyrat.viewer.Widgets import HLine, CropBoxWidget, BoolWidget, FileselWidget, ProductContentWidget_UAVSAR


class Ann(object):
    """Class for loading and interacting with an UAVSAR annotation file."""

    def __init__(self, file):
        """Load in the specified .ann file as a list of strings and
            initialize.

        Arguments:
            file (str): UAVSAR annotation filename to load.

        """
        self.file = file
        fd = open(self.file, 'r')
        self.ann = fd.read().split('\n')
        return

    def query(self, keyword):
        """Query the annotation file for the specified annotation keyword.

        Arguments:
            keyword (str): The keyword to query.

        Returns:
            value: The value of the specified keyword.

        """
        for n in range(len(self.ann)):
            if self.ann[n].startswith(keyword):
                try:
                    val = self.ann[n].rsplit('=')[-1].split(';')[0].split()[0]
                    val = np.array(val, dtype='float')  # if we can convert the string to a number, do so
                    if (val - np.floor(val)) == 0:
                        val = np.array(val,
                                       dtype='int')  # if it's an integer, convert it to one (e.g., number of samples)
                    return val
                except ValueError:  # if we can't convert the string to a number, leave it as a string
                    val = self.ann[n].split('=', maxsplit=1)[-1].split(';')[0].strip()
                    return val

        return None


class UAVSAR(pyrat.ImportWorker):
    """
    Import of UAVSAR SLC product. This class loads the SLC data of one or several bands and / or one
    or several polarisations and / or one of several tracks, together with their meta data into a
    new PyRAT layer(s).

    :param dir: The UAVSAR product directory.
    :type dir: str
    :param bands: Load only this band. '*' to load all bands. Default='*'
    :type band: string
    :param polar: Load only this polarisation. '*' to load all bands. Default='*'
    :type polar: string
    :param tracks: Load only this track. '*' to load all tracks. Default='*'
    :type tracks: string
    :param product: Selects the product component to import. Default='SLC'
    :type product: string
    :param crop: A crop region / subset to import (az_start, az_end, rg_start, rg_end)
    :type crop: tuple
    :param sym: PolSAR symmetrisation. If set, HV and VH are averaged on import. Default=False
    :type sym: bool
    :author:  Gustavo Daniel Martín del Campo Becerra
    """

    gui = {'menu': 'File|Import airborne', 'entry': 'UAVSAR'}
    para = [
        {'var': 'dir', 'value': ''},
        {'var': 'band', 'value': '*'},
        {'var': 'polar', 'value': '*'},
        {'var': 'tracks', 'value': '*'},
        {'var': 'product', 'value': 'SLC'},
        {'var': 'crop', 'value': [0, 0, 0, 0]},
        {'var': 'sym', 'value': False, 'type': 'bool', 'text': 'Cross-polar symmetrisation'}]

    def __init__(self, *args, **kwargs):
        super(UAVSAR, self).__init__(*args, **kwargs)
        self.name = "UAVSAR SLC IMPORT"
        if len(args) == 1:
            self.dir = args[0]

    def reader(self, *args, **kwargs):
        array = []
        meta = []

        code_pos = 5
        code_pos_tracks = 3

        files = glob.glob(os.path.join(self.dir, '*' + '_00' + self.tracks + '_' + '*' + '_' + self.band.upper()
                                       + '*' + self.polar.upper() + '*.slc'))
        files_ann = glob.glob(os.path.join(self.dir, '*' + '_00' + self.tracks + '_' + '*' + '_' + self.band.upper()
                                           + '*' + self.polar.upper() + '*.ann'))
        pols = list(set([os.path.basename(slc).split('_')[code_pos][4:6].upper() for slc in files]))
        bands = list(set([os.path.basename(slc).split('_')[code_pos][0] for slc in files]))
        tracks = list(set([os.path.basename(slc).split('_')[code_pos_tracks][-1].upper() for slc in files]))
        tracks.sort()

        for band in bands:
            if len(bands) > 0 and self.band != '*' and (band not in self.band):
                logging.warning(band + '-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_' + band in f]
                bandfiles_ann = [f for f in files_ann if '_' + band in f]

            trackfiles = []
            trackfiles_ann = []
            for track in tracks:
                trackfiles.append([f for f in bandfiles if '_00' + track in f])
                trackfiles_ann.append([f for f in bandfiles_ann if '_00' + track in f])

            for index in range(len(trackfiles)):
                azsize_crop = self.crop[1] - self.crop[0]
                rgsize_crop = self.crop[3] - self.crop[2]

                if azsize_crop <= 0 or rgsize_crop <= 0:
                    logging.info('')
                    logging.info(bcolors.FAIL + "Incorrect crop size. Aborting." + bcolors.ENDC)
                    return

                barr = np.empty((len(trackfiles[index]), azsize_crop, rgsize_crop), dtype='complex64')
                bmeta = {}
                bmeta['sensor'] = 'UAVSAR'
                bmeta['band'] = band
                bmeta['CH_pol'] = [' '] * len(pols)
                for k, f in enumerate(trackfiles[index]):
                    try:
                        ann = Ann(trackfiles_ann[index][k])
                    except:
                        logging.info('')
                        logging.info(bcolors.FAIL + "Cannot load UAVSAR annotation file. Aborting." + bcolors.ENDC)
                        return

                    bmeta['CH_pol'][k] = os.path.basename(trackfiles[index][k]).split('_')[code_pos][4:6].upper()

                    rngsize_slc = ann.query('slc_1_1x1 Columns')
                    azsize_slc = ann.query('slc_1_1x1 Rows')

                    # Load slc
                    qlimage = np.zeros((azsize_crop, rgsize_crop), dtype='complex64')

                    num_blocks = int(np.sum(azsize_slc) * rngsize_slc * 8 / 1e9)
                    if num_blocks < 2:
                        num_blocks = 2
                    az_vector = np.round(np.linspace(self.crop[0], self.crop[1], num=num_blocks + 1)).astype('int')
                    az_vector[num_blocks] = self.crop[1]

                    for n, azstart in enumerate(az_vector[0:-1]):
                        azend = az_vector[n + 1]
                        azstart_slc = azstart
                        azend_slc = azend
                        seg_start, azoffset_start = self.findsegment(azstart_slc, azsize_slc)
                        seg_end, azoffset_end = self.findsegment(azend_slc, azsize_slc)

                        if seg_start == seg_end:
                            slc = self.getslcblock(f, rngsize_slc, azoffset_start, azoffset_end,
                                                   rngbounds=self.crop[-2:])
                        else:
                            file2 = f[seg_end]
                            slc = self.getslcblock(f, rngsize_slc, azoffset_start, azoffset_end, file2=file2,
                                                   azsize=azsize_slc[seg_start], rngbounds=self.crop[-2:])

                        qlimage[(azstart - self.crop[0]):(azend - self.crop[0]), :] = slc

                    barr[k, ...] = qlimage

                if self.sym is True and barr.ndim == 3 and barr.shape[0] == 4:
                    pol = bmeta['CH_pol']
                    idx_hh = pol.index('HH')
                    idx_vv = pol.index('VV')
                    idx_hv = pol.index('HV')
                    idx_vh = pol.index('VH')
                    barr[idx_hv, ...] = (barr[idx_hv, ...] + barr[idx_vh, ...]) / np.sqrt(2)
                    barr = np.delete(barr, idx_vh, axis=0)
                    bmeta['CH_pol'][idx_hv] = 'XX'
                    bmeta['CH_pol'].remove('VH')
                array.append(barr)
                meta.append(bmeta)

        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], meta[0]
        else:
            return array, meta

    def findsegment(self, az, azsize):
        """For a given azimuth index, return the segment number and
        the azimuth index within the given segment.

        Arguments:
            az: Azimuth index of interest.
            azsize: List containing the azimuth size of each segment.

        Returns:
            seg: Segment number.
            azoff: Azimuth index within the segment.

        """
        azstart = np.insert(np.cumsum(azsize), 0, 0)[0:-1]

        if az <= np.sum(azsize):
            seg = np.max(np.where(azstart <= az))
            azoff = az - azstart[seg]
        else:
            seg = azstart.shape[0] - 1
            azoff = azsize[seg]
            logging.info('')
            logging.info(bcolors.WARNING + 'SLC row index of ' + str(
                az) + ' is larger than the size of the data. Returning maximum index.' + bcolors.ENDC)

        return seg, azoff

    def getslcblock(self, file, rngsize, azstart, azend, rngbounds=None, file2=None,
                    azsize=None):
        """Load SLC data into a NumPy array buffer.  If the file2 argument is
        specified in the arguments, this function will treat the two files as
        consecutive segments, and will join them.

        Arguments:
            file (str): Filename of the first SLC.
            rngsize (int): Number of columns (range bins).  Same for both SLCs.
            azstart (int): Azimuth index at which to start the buffer, in the
                first SLC.
            azend (int): Azimuth index at which to end the buffer.  If file2
                is specified, this is an azimuth index in the second SLC.  If
                file2 is not specified, this is an azimuth index in the first SLC.
                The row specified by azend is not actually included in the buffer,
                as in the Python range() function.  (azend-1) is the last line
                included in the buffer.  To load the entire SLC, azend should be
                equal to the number of rows in SLC.
            rngbounds (tuple, int): Starting and ending range bounds, if range
                subsetting is desired.
            file2 (str): Filename of the second SLC, if combining multiple
                segments is desired.
            azsize (int): Number of rows (azimuth bins) for the first SLC.  Only
                required if file2 is specified.  Otherwise we only load in the
                lines of the SLC before azend.

        Returns:
            block: NumPy array of complex64 datatype, containing the loaded SLC
            data between the specified azimuth bounds.

        """
        if file2 is None:
            byteoffset = rngsize * 8 * azstart
            slc = np.memmap(file, dtype='complex64', mode='c', offset=byteoffset, shape=(azend - azstart, rngsize))
            if rngbounds is not None:
                slc = slc[:, rngbounds[0]:rngbounds[1]]
            return slc
        elif azsize is None:
            logging.info('')
            logging.info(
                bcolors.FAIL + '"file2" argument specified, but "azsize" argument missing. Aborting.' + bcolors.ENDC)
            return
        else:
            byteoffset = rngsize * 8 * azstart
            slca = np.memmap(file, dtype='complex64', mode='c', offset=byteoffset, shape=(azsize - azstart, rngsize))
            slcb = np.memmap(file2, dtype='complex64', mode='c', shape=(azend, rngsize))
            if rngbounds is not None:
                slca = slca[:, rngbounds[0]:rngbounds[1]]
                slcb = slcb[:, rngbounds[0]:rngbounds[1]]
            return np.vstack((slca, slcb))

    @classmethod
    def guirun(cls, viewer):
        para_backup = copy.deepcopy(cls.para)  # keep a deep copy of the default parameters
        wid = UAVSARImportWidget()
        wid.update()
        res = wid.exec_()
        if res == 1:
            plugin = cls(dir=wid.dir, product=wid.product, band=wid.bands, polar=wid.polar, crop=wid.crop, sym=wid.sym,
                         tracks=wid.tracks)
            viewer.statusBar.setMessage(message=plugin.name + ' running', colour='R')
            plugin.run()
            del plugin
            viewer.statusBar.setMessage(message='Ready', colour='G')
            viewer.updateViewer()


class UAVSARImportWidget(QtWidgets.QDialog):
    def __init__(self, parent=None, dir=None):
        super(UAVSARImportWidget, self).__init__(parent)
        self.setWindowTitle("UAVSAR import")
        mainlayout = QtWidgets.QVBoxLayout(self)

        self.dirwidget = FileselWidget(title='UAVSAR product dir', type='opendir')
        self.dirwidget.setvalue(dir)
        mainlayout.addWidget(self.dirwidget)
        mainlayout.addWidget(HLine())
        self.productwidget = ProductContentWidget_UAVSAR(products=['SLC'], tracks=[])
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

        self.bandupdate = lambda: self.update(mode=1)  # update band
        self.productwidget.band.currentIndexChanged.connect(self.bandupdate)

    def update(self, mode=0):
        self.dir = str(self.dirwidget.getvalue())
        self.product = self.productwidget.getvalue(0)
        self.tracks = self.productwidget.getvalue(1)
        self.bands = self.productwidget.getvalue(2)
        self.polar = self.productwidget.getvalue(3)

        code_pos = 5
        code_pos_tracks = 3

        files = glob.glob(os.path.join(self.dir, '*' + '_' + self.bands.upper() + '*' + self.polar.upper() + '*.ann'))

        if mode == 0:
            allfiles = glob.glob(os.path.join(self.dir, '*.ann'))
            self.bands = '*'
            self.polar = '*'
            self.tracks = '*'
        if mode == 1:
            allfiles = glob.glob(os.path.join(self.dir, '*' + '_' + self.bands.upper() + '*.ann'))

        allbands = list(set([os.path.basename(slc).split('_')[code_pos][0] for slc in allfiles]))
        allpols = list(set([os.path.basename(slc).split('_')[code_pos][4:6].upper() for slc in allfiles]))
        alltracks = list(set([os.path.basename(slc).split('_')[code_pos_tracks][-1].upper() for slc in allfiles]))
        alltracks.sort()

        nrg = 0
        naz = 0
        for filename in files:
            try:
                ann = Ann(filename)
            except:
                logging.info('')
                logging.info(bcolors.FAIL + "Cannot load UAVSAR annotation file. Aborting." + bcolors.ENDC)
                return

            # Get SLC dimensions for the 1x1 SLCs
            nrg_max = ann.query('slc_1_1x1 Columns')
            naz_max = ann.query('slc_1_1x1 Rows')
            nrg = max(nrg, nrg_max)
            naz = max(naz, naz_max)
        self.cropwidget.setrange([[0, naz], [0, naz], [0, nrg], [0, nrg]])
        self.cropwidget.setvalues([0, naz, 0, nrg])

        if mode == 0:
            self.productwidget.band.currentIndexChanged.disconnect(self.bandupdate)
            self.productwidget.updatepolar(allpols)
            self.productwidget.updatebands(allbands)
            self.productwidget.updatetracks(alltracks)
            self.productwidget.setvalue(1, self.tracks)
            self.productwidget.setvalue(2, self.bands)
            self.productwidget.setvalue(3, self.polar)
            self.productwidget.band.currentIndexChanged.connect(self.bandupdate)

        elif mode == 1:
            self.productwidget.updatepolar(allpols)
            self.productwidget.setvalue(3, self.polar)

    def accept(self):
        self.product = self.productwidget.getvalue(0)
        self.tracks = self.productwidget.getvalue(1)
        self.bands = self.productwidget.getvalue(2)
        self.polar = self.productwidget.getvalue(3)
        self.crop = self.cropwidget.getvalues()
        self.sym = self.symwidget.getvalue()
        super(UAVSARImportWidget, self).accept()


@pyrat.docstringfrom(UAVSAR)
def uavsar(*args, **kwargs):
    return UAVSAR(*args, **kwargs).run(*args, **kwargs)
