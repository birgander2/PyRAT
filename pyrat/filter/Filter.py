import numpy as np
import scipy as sp
from scipy.ndimage import filters

import pyrat


class Median(pyrat.FilterWorker):
    """
    Median filter of data arrays.

    :author: Joel Amao
    """
    gui = {'menu': 'Tools|Filter', 'entry': 'Median'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']},
        {'var': 'phase', 'value': False, 'type': 'bool', 'text': 'Phase'}
    ]

    def __init__(self, *args, **kwargs):
        super(Median, self).__init__(*args, **kwargs)
        self.name = "MEDIAN FILTER"
        self.blockprocess = True
        self.blockoverlap = self.win[0] // 2 + 1

    def filter(self, array, *args, **kwargs):
        win = self.win
        if array.ndim == 3:
            win = [1] + self.win
        if array.ndim == 4:
            win = [1, 1] + self.win

        if np.iscomplexobj(array):
            return filters.median_filter(array.real, size=win) + 1j * filters.median_filter(array.imag, size=win)
        elif self.phase is True:
            tmp = np.exp(1j * array)
            tmp = filters.uniform_filter(tmp.real, win) + 1j * filters.uniform_filter(tmp.imag, win)
            return np.angle(tmp)
        else:
            return filters.median_filter(array.real, size=win)


@pyrat.docstringfrom(Median)
def median(*args, **kwargs):
    return Median(*args, **kwargs).run(**kwargs)


class Boxcar(pyrat.FilterWorker):
    """
    Boxcar / Moving average (speckle) filter.

    :author: Andreas Reigber
    """

    gui = {'menu': 'Tools|Filter', 'entry': 'Boxcar'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
         'subtext': ['range', 'azimuth']},
        {'var': 'phase', 'value': False, 'type': 'bool', 'text': 'Phase'}
    ]

    def __init__(self, *args, **kwargs):
        super(Boxcar, self).__init__(*args, **kwargs)
        self.name = "BOXCAR FILTER"
        self.blockprocess = True
        self.blockoverlap = self.win[0] // 2 + 1

    def filter(self, array, *args, **kwargs):
        win = self.win
        if array.ndim == 3:
            win = [1] + self.win
        if array.ndim == 4:
            win = [1, 1] + self.win
        array[np.isnan(array)] = 0.0
        if np.iscomplexobj(array):
            return sp.ndimage.filters.uniform_filter(array.real, win) + 1j * sp.ndimage.filters.uniform_filter(
                array.imag, win)
        elif self.phase is True:
            tmp = np.exp(1j * array)
            tmp = sp.ndimage.filters.uniform_filter(tmp.real, win) + 1j * sp.ndimage.filters.uniform_filter(tmp.imag,
                                                                                                            win)
            return np.angle(tmp)
        else:
            return sp.ndimage.filters.uniform_filter(array.real, win)


@pyrat.docstringfrom(Boxcar)
def boxcar(*args, **kwargs):
    return Boxcar(*args, **kwargs).run(**kwargs)


class Gauss(pyrat.FilterWorker):
    """
    Gaussian (speckle) filter...

    :author: Andreas Reigber
    """

    gui = {'menu': 'Tools|Filter', 'entry': 'Gauss'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'float', 'range': [1, 999], 'text': 'Sigma',
         'subtext': ['range', 'azimuth']},
        {'var': 'phase', 'value': False, 'type': 'bool', 'text': 'Phase'}
    ]

    def __init__(self, *args, **kwargs):
        super(Gauss, self).__init__(*args, **kwargs)
        self.name = "GAUSS FILTER"
        self.blockprocess = True
        self.blockoverlap = 2 * int(self.win[0]) + 1

    def filter(self, array, *args, **kwargs):
        win = self.win
        if array.ndim == 3:
            win = [1] + self.win
        if array.ndim == 4:
            win = [1, 1] + self.win
        array[np.isnan(array)] = 0.0
        if np.iscomplexobj(array):
            return sp.ndimage.filters.gaussian_filter(array.real, win) + 1j * sp.ndimage.filters.gaussian_filter(
                array.imag, win)
        elif self.phase is True:
            tmp = np.exp(1j * array)
            tmp = sp.ndimage.filters.gaussian_filter(tmp.real, win) + 1j * sp.ndimage.filters.gaussian_filter(tmp.imag,
                                                                                                              win)
            return np.angle(tmp)
        else:
            return sp.ndimage.filters.gaussian_filter(array.real, win)


@pyrat.docstringfrom(Gauss)
def gauss(*args, **kwargs):
    return Gauss(*args, **kwargs).run(**kwargs)


class PLR(pyrat.FilterWorker):
    """
    Probabilistic Label Relaxation

    To be applied on classification results only (labels in byte / interger representation)

    :author: Andreas Reigber
    """

    gui = {'menu': 'Tools|Filter', 'entry': 'PLR'}
    para = [
        {'var': 'compatibility', 'value': 3, 'type': 'int', 'range': [2, 99], 'text': 'Class compatibility'},
        {'var': 'niter', 'value': 5, 'type': 'int', 'text': 'Number of iterations'}
    ]

    def __init__(self, *args, **kwargs):
        super(PLR, self).__init__(*args, **kwargs)
        self.name = "PLR FILTER"
        self.blockprocess = True
        self.blockoverlap = 3 * self.niter

    def filter(self, array, *args, **kwargs):
        accuracy = 0.9      # initial classification accuracy. Fixed, doesn't influence the result a lot
        win = [3, 3, 1]
        nc = array.max() + 1   # number of classes

        pxx = np.zeros((nc, nc), dtype='f4') + 1.0 / (self.compatibility + nc - 1.0)  # init compatibility matric
        np.fill_diagonal(pxx, self.compatibility / (self.compatibility + nc - 1.0))
        pxx /= pxx.max()

        ms = np.zeros(array.shape + (nc, ), dtype='f4') + (1 - accuracy)   # probability of class membership
        for k in range(nc):
            ms[..., k][array == k] = accuracy                  # initialisation

        for k in range(self.niter):                       # PLR interation
            mss = sp.ndimage.filters.uniform_filter(ms, win)
            q = np.dot(mss, pxx)
            ms *= q
            ms /= np.sum(ms, axis=2)[..., np.newaxis]

        return np.ndarray.astype(np.argmax(ms, axis=2), dtype=array.dtype)

@pyrat.docstringfrom(PLR)
def plr(*args, **kwargs):
    return PLR(*args, **kwargs).run(**kwargs)


class KernelFilter(pyrat.FilterWorker):
    """
    Custom filtering of data arrays given a given (2d) window.

    :author: Joel Amao
    """
    para = [{'var': 'phase', 'value': False, 'type': 'bool', 'text': 'Phase'}]
    def __init__(self, *args, **kwargs):
        super(KernelFilter, self).__init__(*args, **kwargs)
        self.name = "KERNEL FILTER"
        self.blockprocess = True
        self.blockoverlap = 11 // 2 + 1

    def filter(self, array, *args, **kwargs):
        win = self.win
        if array.ndim == 3:
            win = win[np.newaxis,...]
        if array.ndim == 4:
            win = win[np.newaxis, np.newaxis, ...]

        if np.iscomplexobj(array):
            return filters.convolve(array.real, weights=win) + 1j * filters.convolve(array.imag, weights=win)
        elif self.phase is True:
            tmp = np.exp(1j * array)
            tmp = filters.convolve(tmp.real, weights=win) + 1j * filters.convolve(tmp.imag, weights=win)
            return np.angle(tmp)
        else:
            return filters.convolve(array.real, weights=win)

@pyrat.docstringfrom(KernelFilter)
def kernelfilter(*args, **kwargs):
    return KernelFilter(*args, **kwargs).run(**kwargs)

