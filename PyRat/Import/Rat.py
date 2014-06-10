import PyRat
import STEtools as STE
import numpy as np
import pdb

class Rat(PyRat.ImportWorker):
    def __init__(self, *args, **kwargs):
        super(Rat, self).__init__(*args, **kwargs)    
        self.name = "RAT IMPORT"
        
    def reader(self, *args, **kwargs):
        array = STE.rrat(self.filename)
        
# pixel interleaved vector -> band interleaved vector        
        
        if array.ndim == 3 and array.shape[2] < array.shape[0] and array.shape[2] < array.shape[1]:
            array = np.rollaxis(array,2)
        
# pixel interleaved matrix -> band interleaved vector        
        
        if array.ndim == 4 and array.shape[2] < array.shape[0] and array.shape[2] < array.shape[1] and array.shape[3] < array.shape[0] and array.shape[3] < array.shape[1]:
            array = np.rollaxis(np.rollaxis(array,2),3,start=1)
            array = np.reshape(array,(array.shape[0]*array.shape[1],array.shape[2],array.shape[3]))
        
        return array, None
