import PyRat
import numpy as np

import ipdb; stop = ipdb.set_trace

class Rotate(PyRat.FilterWorker):
    """
    Rotate data set by multiples of 90 degrees
    
    :param times: How often to rotate
    :type int:
    :author: Andreas Reigber
    """
    def __init__(self, *args, **kwargs):
        super(Rotate, self).__init__(*args, **kwargs)    
        self.name = "ROTATE"
        self.times = 1
        
    def filter(self, array, *args, **kwargs):
        if array.ndim == 3: array = np.rollaxis(array, axis=0, start=3)
        array = np.rot90(array, self.times)
        if array.ndim == 3: array = np.rollaxis(array, axis=2)
        return array
