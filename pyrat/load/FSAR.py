import pyrat
import glob, os, logging
import numpy as np
import pdb

class FSAR_SLC(pyrat.ImportWorker):
    """
    Import of DLR F-SAR SLC product
    
    :author: Andreas Reigber
    """
    def __init__(self, *args, **kwargs):
        super(FSAR_SLC, self).__init__(*args, **kwargs)    
        self.name = "FSAR SLC IMPORT"
        
    def reader(self, *args, **kwargs):
        meta  = {}
        meta['sensor'] = 'DLR F-SAR'
        array = False
        
        files = glob.glob(self.filename+'/RGI/RGI-SR/slc*rat')
        
        pdb.set_trace()
        return array, meta
