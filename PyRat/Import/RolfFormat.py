import PyRat
import numpy as np

import ipdb; stop = ipdb.set_trace

class RolfFormat(PyRat.ImportWorker):
    """
    Import of Rolf's old binary format with 2 longs as header
    
    :param filename: The filename of the DCSLC file to import.
    :type array: string
    :param dtype: The dtype of the array to import. Default: 'float32'
    :type array: string
    :author: Andreas Reigber
   
    """
    def __init__(self, *args, **kwargs):
        super(RolfFormat, self).__init__(*args, **kwargs)    
        self.name = "OLD ROLF FORMAT IMPORT"
        self.dtype = 'float32'
        
    def reader(self, *args, **kwargs):
        
        lun = open(self.filename,'rb')
        nrx, nry = np.fromfile(lun, dtype=np.uint32, count=2).byteswap()
        array = np.empty((nry,nrx),dtype=self.dtype)
        lun.readinto(array.data)
        array = array.byteswap()
        lun.close()
        
        return array, None
