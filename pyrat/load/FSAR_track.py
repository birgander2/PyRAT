import pyrat, glob, os, logging
import numpy as np
from pyrat.load.tools import RatFile, Xml2Py


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
        self.name  = "FSAR REAL/REF TRACK IMPORT"
        self.block = 'T'
        if 'crop' not in self.__dict__:  self.crop = (0,0,0,0)

    def reader(self, *args, **kwargs):
        if self.product == 'RGI-SLC':
            #head = 'reftr_sar_resa'
            head = '_sar_resa'
            src = ('RGI','RGI-TRACK')
        if self.product == 'RGI-AMP':
            #head = 'reftr_sar_resa'
            head = '_sar_resa'
            src = ('RGI','RGI-TRACK')
        if self.product == 'INF-SLC':
            #head = 'reftr_sar_resa'
            head = '_sar_resa'
            src = ('INF','INF-TRACK')

        files = glob.glob(os.path.join(self.dir, src[0], src[1], '*'+head+'*'+self.bands.upper()+self.polarisations.lower()+'*.rat'))
        bands = list(set([os.path.basename(slc).split('_')[4][0] for slc in files]))
        pols = list(set([os.path.basename(slc).split('_')[4][1:3] for slc in files]))
        tracks = list(set([os.path.basename(slc).split('_')[0][0:6] for slc in files])) # contains the path to the 2 track files (reference and real)

        array = []
        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band+'-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_'+band in f] # list of files from the same band
                fil = RatFile(bandfiles[0])
                naz = fil.dim[1]
                block = list(self.crop)
                block[2] = 0
                block[3] = 7
                if block[1] == 0:
                    block[1] = naz
                daz = block[1]-block[0]
                drg = block[3]-block[2]

                barr = np.empty((len(bandfiles)/2, 7, daz))

                for k, pol in enumerate(pols):
                    polfiles = [f for f in bandfiles if pol+'_' in f]
                    for i, f in enumerate(polfiles):
                        logging.info("Found "+f)
                        if 'reftr' in f:
                            barr[k, 0:3, ...] = RatFile(f).read(block=(2*block[0], 1, 2*daz, 3), step=2)
                            barr[k, 6, ...] = RatFile(f).read(block=(2*block[0], 0, 2*daz, 1), step=2)
                        elif 'track' in f:
                            barr[k, 3:6, ...] = RatFile(f).read(block=(2*block[0], 1, 2*daz, 3), step=2)
                array.append(barr)
        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0].T, None
        else:
            return array, None

'''
        # TODO: use PyRAT's rrat instead of the one from STEtools
        #track = RatFile(file[0])
        track = STEtools.rrat(file[0])
        head = 'reftr_sar_resa'
        file = glob.glob(os.path.join(self.dir, src[0], src[1], head+'*'+self.bands.upper()+self.polarisations.lower()+'*.rat'))
        reftr = STEtools.rrat(file[0])

        meta_track = np.empty([7, daz])
        meta_track[0:3,:] = reftr[1:,2*block[0]:2*block[1]:2]
        meta_track[3:6,:] = track[1:,2*block[0]:2*block[1]:2]
        meta_track[6,:] = reftr[0:1,2*block[0]:2*block[1]:2]    # load time
        return np.transpose(meta_track), None
'''

def fsar_track(*args, **kwargs):
    return FSAR_track(*args, **kwargs).run(**kwargs)