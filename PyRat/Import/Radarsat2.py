import PyRat
import numpy as np
from osgeo import gdal
import pdb

class Radarsat2(PyRat.ImportWorker):
    def __init__(self, *args, **kwargs):
        super(Radarsat2, self).__init__(*args, **kwargs)    
        self.name = "RADARSAT2 IMPORT"
        
    def reader(self, *args, **kwargs):
        attrs = kwargs['attrs']
        
        ds = gdal.Open(self.filename)
        array = []
        meta  = ds.GetMetadata()
        chn_attrs = []
        for band in range(ds.RasterCount):
            Band = ds.GetRasterBand(band+1)
            array.append(np.array(Band.ReadAsArray()))
            chn_attrs.append(Band.GetMetadata())
        ds = None
        n_channels = len(array)
        array = np.rollaxis(np.dstack(array),2)
        attrs['CH_pol'] = []
        for k in range(n_channels):
            attrs['CH_pol'].append(chn_attrs[k]['POLARIMETRIC_INTERP'])  
        for k,v in meta.items():
            attrs[k] = v
        return array  
