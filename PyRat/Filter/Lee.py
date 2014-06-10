import PyRat
import scipy,pdb

class Lee(PyRat.FilterWorker):
    """
    Lee's classical speckle filter from 1981. Not the best one...

    :author: Andreas Reigber
    :param array: The image to filter (2D numpy.ndarray)
    :type array: float
    :param box: The filter window size
    :type arr2: integer
    :param looks=1.0: The effective number of looks of the input image.
    :type looks: float
    :returns: filtered image
    """

    def __init__(self, *args, **kwargs):
        super(Lee, self).__init__(*args, **kwargs)    
        self.name = "LEE FILTER"
        self.blockprocess  = True
        self.allowed_dtype = ['float32','float64']
        self.allowed_ndim  = [2]
        if 'win'   not in self.__dict__: self.win=7
        if 'looks' not in self.__dict__: self.looks=1
       
        
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

