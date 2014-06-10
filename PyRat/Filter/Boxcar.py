import PyRat
import scipy, numpy
from scipy.ndimage import filters
import pdb, logging, time

class Boxcar(PyRat.FilterWorker):
    """
    Simple Boxcar filter...
    """

    def __init__(self, *args, **kwargs):
        super(Boxcar, self).__init__(*args, **kwargs)    
        self.name = "BOXCAR FILTER"

        if 'phase' not in self.__dict__: self.phase = False
        if 'win'   not in self.__dict__: 
            self.win=[7,7]
        elif isinstance(self.win, int): 
            self.win = [self.win]*2
        
        self.blockprocess = True
        self.blockoverlap = self.win[0]/2+1
        
    def filter(self, array, *args, **kwargs):
        win = self.win
        if array.ndim == 3: win = [1] + self.win
        if array.ndim == 4: win = [1, 1] + self.win
        if numpy.iscomplexobj(array):
            return(filters.uniform_filter(array.real,win) + 1j * filters.uniform_filter(array.imag,win))
        elif self.phase == True:
            return(numpy.angle(smooth(numpy.exp(1j*array),win)))
        else:
            return(filters.uniform_filter(array.real,win))

