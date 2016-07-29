import pyrat
import numpy as np
from osgeo import gdal
import logging


class GenericGdal(pyrat.ImportWorker):
    """
    Generic GDAL format reader.
    Only main meta information is imported, not the one of the individual bands.

    :author: Andreas Reigber
    """

    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'Product file'}]

    def __init__(self, *args, **kwargs):
        super(GenericGdal, self).__init__(*args, **kwargs)
        self.name = "GENERIC GDAL IMPORT"
        self.vblock = True
        if len(args) == 1:
            self.file = args[0]

    def getsize(self, *args, **kwargs):
        self.ds = gdal.Open(self.file)
        if self.ds is not None:
            self.band = []
            for band in range(self.ds.RasterCount):
                self.band.append(self.ds.GetRasterBand(band + 1))
            return self.ds.RasterCount, self.ds.RasterXSize, self.ds.RasterYSize
        else:
            logging.error("ERROR: GDAL format not recognised!")
            return False, False

    def block_reader(self, *args, **kwargs):
        array = []
        for band in self.band: array.append(band.ReadAsArray(xoff=0, yoff=kwargs['block'][2], win_ysize=self.blocksize))
        out = np.empty((len(array),) + array[0].shape[::-1], dtype=array[0].dtype)
        for k in range(len(array)):
            out[k, ...] = np.rot90(array[k])
        return out

    def close(self, *args, **kwargs):
        self.ds = None  # correct according to GDAL manual!!??

    def getmeta(self, *args, **kwargs):
        meta = {}
        metain = self.ds.GetMetadata()
        meta.update(metain)
        return meta


@pyrat.docstringfrom(GenericGdal)
def genericgdal(*args, **kwargs):
    return GenericGdal(*args, **kwargs).run(*args, **kwargs)


class GeoTiff(GenericGdal):
    """
    GeoTiff format reader (using GDAL)
    """
    def __init__(self, *args, **kwargs):
        super(GeoTiff, self).__init__(*args, **kwargs)
        self.name = "GEOTIFF IMPORT"


@pyrat.docstringfrom(GeoTiff)
def geotiff(*args, **kwargs):
    return GeoTiff(*args, **kwargs).run(*args, **kwargs)

