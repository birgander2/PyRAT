from __future__ import print_function
import pyrat
import numpy as np
from skimage import filter

class Canny(pyrat.FilterWorker):
    """
    Canny edge detector

    :author: Andreas Reigber
    :param sigma: The standard deviation of the filter (higher: less sensitive to noise)
    :type sigma: float

    """

    gui = {'menu': 'SAR|Edge detection', 'entry': 'Canny'}
    para = [{'var': 'sigma', 'value': 2.0, 'type': 'float',  'range': [0.1, 999.], 'text': 'Sigma'}]

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.blockprocess = True
        self.blockoverlap = int(self.sigma)*4

    def filter(self, array, *args, **kwargs):
        array[np.isnan(array)] = 0.0
        array = filter.canny(array, sigma=self.sigma)
        return array


def canny(*args, **kwargs):
    return Canny(*args, **kwargs).run(**kwargs)
