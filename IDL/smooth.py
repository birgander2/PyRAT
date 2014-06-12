import numpy
from scipy.ndimage import filters
def smooth(array, box, phase=False):
    """
    Imitates IDL's smooth function. Can also (correctly) smooth interferometric phases with the phase=True keyword.
    """
    if numpy.iscomplexobj(array):
        return filters.uniform_filter(array.real,box) + 1j * filters.uniform_filter(array.imag,box)
    elif phase == True:
        return numpy.angle(smooth(numpy.exp(1j*array),box))
    else:
       return filters.uniform_filter(array.real,box)

