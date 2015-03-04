import pyrat, os
import numpy as np
from osgeo import gdal


class TSX(pyrat.ImportWorker):
    """
    Import of TSX/TDX satellite data
    """

    gui = {'menu': 'File|Open external', 'entry': 'TerraSAR-X'}
    para = [{'var': 'filename', 'value': '', 'type': 'openfile', 'text': 'Product xml file'}]

    def __init__(self, *args, **kwargs):
        super(TSX, self).__init__(*args, **kwargs)
        self.name = "TSX / TDX IMPORT"
        
    def reader(self, *args, **kwargs):
        meta  = {}
        meta['sensor'] = 'TSX'
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
            array = [np.squeeze(np.rollaxis(np.dstack(array),2))]
            
            meta['CH_pol'] = []
            for k in range(n_channels):
                meta['CH_pol'].append(chn_attrs[k]['POLARIMETRIC_INTERP'])
        else:
            logging.error("File not found: "+self.filename)
        return array , meta
