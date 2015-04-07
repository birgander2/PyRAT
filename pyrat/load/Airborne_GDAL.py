import pyrat, os
import logging
from osgeo import gdal
import numpy as np

class AIRSAR(pyrat.ImportWorker):
    """
    Import of JPL AIRSAR covariance files (\*.STK).

    **author:** Andreas Reigber\n
    **status:** --beta-- Mostly untested.
    """

    gui = {'menu': 'File|Import airborne', 'entry': 'AIRSAR'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'Product file'}]

    def __init__(self, *args, **kwargs):
        super(AIRSAR, self).__init__(*args, **kwargs)
        self.name = "AIRSAR IMPORT"

    def getsize(self, *args, **kwargs):
        self.ds = gdal.Open(self.file)
        if self.ds is not None:
            self.band = []
            for band in range(self.ds.RasterCount):
                self.band.append(self.ds.GetRasterBand(band+1))
            return 3, 3, self.ds.RasterYSize, self.ds.RasterXSize
        else:
            logging.error("ERROR: product not recognised!")
            return False, False

    def block_reader(self, *args, **kwargs):
        array = []

        for band in self.band:
            array.append(band.ReadAsArray(xoff=0, yoff=kwargs['block'][0], win_ysize=self.blocksize))

        out = np.empty((3, 3)+array[0].shape, dtype=array[0].dtype)
        out[0, 0, ...] = array[0] + 0j   # HH-HH
        out[1, 1, ...] = array[5] + 0j   # VV-VV
        out[2, 2, ...] = array[3] + 0j   # 2* XX-XX
        out[0, 1, ...] = array[2]   # sqrt(2) HH-VV
        out[1, 0, ...] = np.conj(array[2])
        out[0, 2, ...] = array[1]   # sqrt(2) HH-XX
        out[2, 0, ...] = np.conj(array[1])
        out[1, 2, ...] = np.conj(array[4])
        out[2, 1, ...] = array[4]   # sqrt(2) HV-VV

        return out

    def close(self, *args, **kwargs):
        self.ds = None          # correct according to GDAL manual!!??

    def getmeta(self, *args, **kwargs):
        meta = {}
        metain = self.ds.GetMetadata()
        meta.update(metain)
        meta['CH_pol'] = ['HHHH*', 'HHVV*', 'HHXX*', 'VVHH*', 'VVVV*', 'VVXX*', 'XXHH*', 'XXVV*', 'XXXX*']
        return meta


class Convair(pyrat.ImportWorker):
    """
    Import of CONVAIR-580 data

    **author:** Andreas Reigber\n
    **status:** --beta-- Mostly untested.
    """

    gui = {'menu': 'File|Import airborne', 'entry': 'CONVAIR'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'Product file'}]

    def __init__(self, *args, **kwargs):
        super(Convair, self).__init__(*args, **kwargs)
        self.name = "CONVAIR IMPORT"
        self.vblock = True

    def getsize(self, *args, **kwargs):
        self.ds = gdal.Open(self.file)
        if self.ds is not None:
            self.band = []
            for band in range(self.ds.RasterCount):
                self.band.append(self.ds.GetRasterBand(band+1))
            return self.ds.RasterCount, self.ds.RasterXSize, self.ds.RasterYSize
        else:
            logging.error("ERROR: product not recognised!")
            return False, False

    def block_reader(self, *args, **kwargs):
        array = []
        for band in self.band: array.append(band.ReadAsArray(xoff=0, yoff=kwargs['block'][2], win_ysize=self.blocksize))
        out = np.empty((len(array), )+array[0].shape[::-1], dtype=array[0].dtype)
        for k in range(len(array)):
            out[k, ...] = np.rot90(array[k])
        return out

    def close(self, *args, **kwargs):
        self.ds = None          # correct according to GDAL manual!!??

    def getmeta(self, *args, **kwargs):
        meta = {}
        metain = self.ds.GetMetadata()
        meta.update(metain)

        meta['CH_pol'] = []
        for band in self.band:
            metain = band.GetMetadata()
            meta['CH_pol'].append(metain['POLARIMETRIC_INTERP'].upper())
        return meta
