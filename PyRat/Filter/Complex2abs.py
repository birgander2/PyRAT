import PyRat
import numpy as np

class Complex2abs(PyRat.FilterWorker):
    """
    Complex to absolute conversion
    """

    def __init__(self, *args, **kwargs):
        super(Complex2abs, self).__init__(*args, **kwargs)    
        self.name = "Complex -> Abs"
        self.blockprocess = True
        
    def filter(self, array, *args, **kwargs):
        return np.abs(array)
  
