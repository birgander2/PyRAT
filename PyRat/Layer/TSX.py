import logging
import PyRat
import numpy as np
from osgeo import gdal
from osgeo.gdalconst import *
import pdb

class TSX(PyRat.LayerWorker):
    def __init__(self, filename=None, *args, **kwargs):
        self.ds = gdal.Open(filename, GA_ReadOnly)
        if self.ds is None:
            logging.error('Could not open file : '+filename)
            return
        super(TSX, self).__init__(*args, **kwargs)    
        nx = self.ds.RasterXSize
        ny = self.ds.RasterYSize
        nb = self.ds.RasterCount
        
        self.attrs = {}
        self.attrs['_type']   = 'TSX/TDX'
        self.attrs['_shape']  = (nb, ny, nx)
        self.attrs['_lshape'] = (nb,)
        self.attrs['_dshape'] = (ny, nx)
        self.attrs['_dtype']  = None
                
    def getMeta(self, key=False):
        meta = self.ds.GetMetadata()
        return meta
    
    def getData(self, block=[0,0,0,0], layer=None):
        if block[1] == 0: block[1] =  self.attrs['_dshape'][0]
        if block[3] == 0: block[3] =  self.attrs['_dshape'][1]
        if layer == None or layer == self.name:
            array = []
            for band in range(self.ds.RasterCount):
                array.append(np.array(self.ds.GetRasterBand(band+1).ReadAsArray()))
            return np.squeeze(np.rollaxis(np.dstack(array),2))
        elif 'D' in layer:
            band = self.ds.GetRasterBand(int(layer.split('/')[2][1:]))
            return np.squeeze(band.ReadAsArray(block[0],block[2],block[1]-block[0],block[3]-block[2]))
        else:
            logging.error('Layer name unknown') 
            pdb.set_trace()

