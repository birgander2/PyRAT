import PyRat, os
import numpy as np
from osgeo import gdal
import pdb

class Radarsat2(PyRat.ImportWorker):
    def __init__(self, *args, **kwargs):
        super(Radarsat2, self).__init__(*args, **kwargs)    
        self.name = "RADARSAT2 IMPORT"
        
    def reader(self, *args, **kwargs):
        meta  = {}
        meta['sensor'] = 'Radarsat 2'
        array = None
        
        if os.path.isfile(self.filename):
            ds = gdal.Open(self.filename)
            array = []
            metain  = ds.GetMetadata()
            meta.update(metain)
            
            chn_attrs = []
            for band in range(ds.RasterCount):
                Band = ds.GetRasterBand(band+1)
                array.append(np.array(Band.ReadAsArray()))
                chn_attrs.append(Band.GetMetadata())
            ds = None
            n_channels = len(array)
            array = np.rollaxis(np.dstack(array),2)
            
            meta['CH_pol'] = []
            for k in range(n_channels):
                meta['CH_pol'].append(chn_attrs[k]['POLARIMETRIC_INTERP'])  
        else:
            logging.error("File not found: "+self.filename)

        return array , meta
