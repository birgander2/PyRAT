import scipy
import numpy

def lee(array, win=7, looks=1.0):
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
    sig2  = 1.0 / looks
    sfak  = 1.0 + sig2
    m2arr  = scipy.ndimage.filters.uniform_filter(array**2, size=win)
    marr   = scipy.ndimage.filters.uniform_filter(array, size=win)
    vary   = (m2arr - marr**2).clip(1e-10)
    varx   = ((vary - marr**2*sig2)/sfak).clip(0)
    k      = varx / vary
    out    = marr + (array-marr) * k
    return out

