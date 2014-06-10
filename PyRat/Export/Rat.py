import PyRat
import STEtools as STE
import numpy as np
import pdb

class Rat(PyRat.ExportWorker):
    def __init__(self, *args, **kwargs):
        super(Rat, self).__init__(*args, **kwargs)    
        self.name = "RAT EXPORT"
       
    def writer(self, array, *args, **kwargs):
        STE.srat(self.filename, np.squeeze(array))
        return True
