import pyrat
import glob, os
import logging
import copy
import numpy as np
from PyQt4 import QtGui, QtCore
from pyrat.load.tools import RatFile, Xml2Py
from pyrat.viewer.Dialogs import FlexFilesel
from pyrat.viewer.Widgets import HLine, CropBoxWidget, FileselWidget, ProductContentWidget


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
        if 'crop' not in self.__dict__:
            self.crop = (0, 0, 0, 0)

    def reader(self, *args, **kwargs):

        head = 'slantdem_full'
        src = ('RGI','RGI-AUX')

        files = glob.glob(os.path.join(self.dir, src[0], src[1], head+'*'+self.bands.upper()+'*.rat'))
        bands = list(set([os.path.basename(slc).split('_')[2][0] for slc in files]))
        pols = list(set([os.path.basename(slc).split('_')[2][1:3] for slc in files]))

        # stop()
        # match = self.match if hasattr(self, 'match') else ''
        # files = glob.glob(self.filename+'/RGI/RGI-SR/slc_*'+match+'_*rat')
        # bands = set([f.split('_')[-2][0] for f in files])

        array = []
        meta = []
        meta.append({})
        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band+'-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_'+band in f]
                fil = RatFile(bandfiles[0])
                naz = fil.dim[1]
                nrg = fil.dim[0]
                block = list(self.crop)
                if block[1] == 0:
                    block[1] = naz
                if block[3] == 0:
                    block[3] = nrg
                daz = block[1]-block[0]
                drg = block[3]-block[2]

                barr = np.empty((len(bandfiles), daz, drg), dtype='complex64')
                for k, f in enumerate(bandfiles):
                    logging.info("Found "+f)
                    barr[k, ...] = RatFile(f).read(block=(block[2], block[0], drg, daz))
                array.append(barr)

        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], meta[0]
        else:
            return array, meta

def fsar_dem(*args, **kwargs):
    return FSAR_dem(*args, **kwargs).run(**kwargs)

