import PyRat
import scipy,numpy
from scipy.ndimage import filters
import pdb, logging, time

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
class Lee(PyRat.FilterWorker):
    """
    Lee's classical speckle filter from 1981. Not the best one...

    :author: Andreas Reigber
    """

    def __init__(self, *args, **kwargs):
        super(Lee, self).__init__(*args, **kwargs)    
        self.name = "LEE FILTER"
        self.allowed_dtype = ['float32','float64']
        self.allowed_ndim  = [2]
        if 'win'   not in self.__dict__: self.win=7
        if 'looks' not in self.__dict__: self.looks=1
        self.blockprocess = True
        self.blockoverlap = self.win[0]/2+1
       
        
    def filter(self, array, *args, **kwargs):
        sig2  = 1.0 / self.looks
        sfak  = 1.0 + sig2
        m2arr  = scipy.ndimage.filters.uniform_filter(array**2, size=self.win)
        marr   = scipy.ndimage.filters.uniform_filter(array, size=self.win)
        vary   = (m2arr - marr**2).clip(1e-10)
        varx   = ((vary - marr**2*sig2)/sfak).clip(0)
        k      = varx / vary
        out    = marr + (array-marr) * k
        return out

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
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
            return filters.uniform_filter(array.real,win) + 1j * filters.uniform_filter(array.imag,win)
        elif self.phase == True:
            tmp = numpy.exp(1j*array)
            tmp = filters.uniform_filter(tmp.real,win) + 1j * filters.uniform_filter(tmp.imag,win)
            return numpy.angle(tmp)
        else:
            return filters.uniform_filter(array.real,win)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
class Gauss(PyRat.FilterWorker):
    """
    Simple Gaussian filter...
    """

    def __init__(self, *args, **kwargs):
        super(Gauss, self).__init__(*args, **kwargs)    
        self.name = "GAUSS FILTER"

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
            return filters.gaussian_filter(array.real,win) + 1j * filters.gaussian_filter(array.imag,win)
        elif self.phase == True:
            tmp = numpy.exp(1j*array)
            tmp = filters.gaussian_filter(tmp.real,win) + 1j * filters.gaussian_filter(tmp.imag,win)
            return numpy.angle(tmp)
        else:
            return filters.gaussian_filter(array.real,win)

