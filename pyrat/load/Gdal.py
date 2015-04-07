import pyrat
import numpy as np
from osgeo import gdal


class Gdal(pyrat.ImportWorker):
    def __init__(self, *args, **kwargs):
        super(Gdal, self).__init__(*args, **kwargs)
        self.name = "GDAL IMPORT"
        
    def reader(self, *args, **kwargs):
        ds = gdal.Open(self.filename)
        array = []
        meta = ds.GetMetadata()
        for band in range(ds.RasterCount):
            array.append(np.array(ds.GetRasterBand(band+1).ReadAsArray()))
        ds = None
        return array, meta  # Workaround for missing multichannel support
