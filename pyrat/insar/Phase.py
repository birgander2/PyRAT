import pyrat
import numpy as np


class Phase(pyrat.FilterWorker):
    """
    Calc InSAR phase between two images
    """

    def __init__(self, *args, **kwargs):
        super(Phase, self).__init__(*args, **kwargs)    
        self.name = "CALC InSAR PHASE"
        self.blockprocess = True
        self.nthreads = 1
        
    def filter(self, array, *args, **kwargs):
        return np.angle(array[0] * np.conj(array[1]))
