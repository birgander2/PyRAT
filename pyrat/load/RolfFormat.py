import pyrat
import numpy as np


class RolfFormat(pyrat.ImportWorker):
    """
    Import of Rolf's old binary format with 2 longs as header
    
    :param filename: The filename of the DCSLC file to import.
    :type array: string
    :param crop: Sets a crop area to read. Default: (0,0,0,0) = entire file
    :type array: tuple
    :param dtype: The dtype of the array to import. Default: 'float32'
    :type array: string
    :author: Andreas Reigber
   
    """
    def __init__(self, *args, **kwargs):
        super(RolfFormat, self).__init__(*args, **kwargs)    
        self.name = "OLD ROLF FORMAT IMPORT"
        if 'dtype' not in self.__dict__: self.dtype = 'f4'
        if 'crop' not in self.__dict__:  self.crop = (0,0,0,0)
        
    def reader(self, *args, **kwargs):
        
        lun = open(self.filename,'rb')
        nrx, nry = np.fromfile(lun, dtype=np.uint32, count=2).byteswap()
        block  = list(self.crop)
        if block[1] == 0: block[1] = nry
        if block[3] == 0: block[3] = nrx
        
        offset = block[0] * nrx * np.dtype(self.dtype).itemsize
        lun.seek(offset, 1)
        
        dy = block[1]-block[0]
        array = np.empty((dy,nrx),dtype=self.dtype)
        lun.readinto(array.data)
        array = array[:,block[2]:block[3]]
        array = array.byteswap()
        lun.close()
        
        return array, None

def rolf(*args, **kwargs):
    return RolfFormat(*args, **kwargs).run(**kwargs)
