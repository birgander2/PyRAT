import pyrat, os
import numpy as np
import logging
from osgeo import gdal


class Radarsat2(pyrat.ImportWorker):
    gui = {'menu': 'File|Open external', 'entry': 'RadarSAR-2'}
    para = [{'var': 'filename', 'value': '', 'type': 'openfile', 'text': 'Product xml file'}]

    def __init__(self, *args, **kwargs):
        super(Radarsat2, self).__init__(*args, **kwargs)    
        self.name = "RADARSAT2 IMPORT"
        
    def reader(self, *args, **kwargs):
        meta  = {}
        meta['sensor'] = 'Radarsat-2'
        array = None
        
        if os.path.isfile(self.filename):
            ds = gdal.Open(self.filename)
            array = []
            metain = ds.GetMetadata()
            meta.update(metain)
            
            chn_attrs = []
            for band in range(ds.RasterCount):
                band = ds.GetRasterBand(band+1)
                array.append(np.array(band.ReadAsArray()))
                chn_attrs.append(band.GetMetadata())
            ds = None
            n_channels = len(array)
            array = np.rollaxis(np.dstack(array),2)
            
            meta['CH_pol'] = []
            for k in range(n_channels):
                meta['CH_pol'].append(chn_attrs[k]['POLARIMETRIC_INTERP'])  
        else:
            logging.error("File not found: "+self.filename)

        return array, meta
