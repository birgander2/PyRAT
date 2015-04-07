from __future__ import print_function
import pyrat
import numpy as np
from scipy import ndimage
import logging


class CoVar(pyrat.FilterWorker):
    """
    Coefficient of variation calculation

    This is a simple SAR image texture descriptor. It should be applied on intensity images!
    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Texture', 'entry': 'Coefficient of variation'}
    para = [{'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']}]

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.blockprocess = True
        self.blockoverlap = self.win[0]//2 + 1

    def filter(self, array, *args, **kwargs):
        shp = array.shape
        array = array.reshape((np.prod(shp[0:-2]),) + shp[-2:])   # reshape to (nch, ny, nx)

        for k, amp in enumerate(np.abs(array)):  # loop over channels
            array[k, ...] = self.covar(amp, self.win)

        array = array.reshape(shp)      # back to original shape
        return array

    @staticmethod
    def covar(amp, win):
        np.seterr(divide='ignore', invalid='ignore')
        samp = ndimage.filters.uniform_filter(amp, size=win)**2
        samp = (ndimage.filters.uniform_filter(amp**2, size=win) - samp)/samp
        samp[~np.isfinite(samp)] = 0.0
        np.seterr(divide='warn', invalid='warn')
        return samp


def covar(*args, **kwargs):
    return CoVar(*args, **kwargs).run(**kwargs)

############################################################################################
############################################################################################
############################################################################################


class Inhomogenity(pyrat.FilterWorker):
    """
    Texture inhomogenenity

    This is a simple SAR image texture descriptor. It should be applied on intensity images!
    The function estimates the k-factor of the Lee speckle filter and uses it as a descriptor
    for local image inhomogenity. Amplitude damping: See P. Leducq, PhD, University of Rennes, 2006
    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Texture', 'entry': 'Inhomogenity'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']},
        {'var': 'looks', 'value': 1.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
        {'var': 'damp', 'value': False, 'type': 'bool', 'text': 'amplitude damping'}
        ]

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.blockprocess = True
        self.blockoverlap = self.win[0]//2 + 1

    def filter(self, array, *args, **kwargs):
        shp = array.shape
        array = array.reshape((np.prod(shp[0:-2]),) + shp[-2:])   # reshape to (nch, ny, nx)

        for k, amp in enumerate(np.abs(array)):  # loop over channels
            array[k, ...] = self.inhomogenity(amp, self.looks, self.win, damp=self.damp)

        array = array.reshape(shp)      # back to original shape
        return array

    @staticmethod
    def inhomogenity(amp, looks, win, damp=False):
        sig2 = 1.0 / looks
        sfak = 1.0 + sig2
        m2arr = ndimage.filters.uniform_filter(amp**2, size=win)
        marr = ndimage.filters.uniform_filter(amp, size=win)
        vary = (m2arr - marr ** 2).clip(1e-10)
        varx = ((vary - marr ** 2 * sig2) / sfak).clip(0)
        kfac = varx / vary
        if damp is True:
            kfac *= ((amp-marr)/(amp+marr)).clip(0)
        return kfac


def inhomgenenity(*args, **kwargs):
    return Inhomogenity(*args, **kwargs).run(**kwargs)

############################################################################################
############################################################################################
############################################################################################

try:
    from skimage.filter.rank import entropy as skentropy
    from pyrat.viewer.tools import sarscale

    class Entropy(pyrat.FilterWorker):
        """
        Local image entropy

       :author: Andreas Reigber
        """

        gui = {'menu': 'SAR|Texture', 'entry': 'Entropy'}
        para = [{'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']}]

        def __init__(self, *args, **kwargs):
            pyrat.FilterWorker.__init__(self, *args, **kwargs)
            self.blockprocess = True
            self.blockoverlap = self.win[0]//2 + 1

        def filter(self, array, *args, **kwargs):
            shp = array.shape
            array = array.reshape((np.prod(shp[0:-2]),) + shp[-2:])   # reshape to (nch, ny, nx)
            selem = np.ones(self.win)
            for k, amp in enumerate(np.abs(array)):  # loop over channels
                array[k, ...] = skentropy(sarscale(amp), selem)

            array = array.reshape(shp)      # back to original shape
            return array

    def entropy(*args, **kwargs):
        return Entropy(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("Entropy module requires skimage library (install?)")

########################################################################################################################
########################################################################################################################
########################################################################################################################

try:
    from skimage.feature import greycomatrix, greycoprops
    from pyrat.viewer.tools import sarscale
    from numpy.lib.stride_tricks import as_strided

    class GLCM(pyrat.FilterWorker):
        """
        Grey level co-occurence matrix features *BETA*

       :author: Andreas Reigber
        """

        gui = {'menu': 'SAR|Texture', 'entry': 'GLCM'}
        para = [
            {'var': 'win', 'value': 7, 'type': 'int', 'range': [3, 999], 'text': 'Window size'},
            {'var': 'levels', 'value': 16, 'type': 'int', 'range': [8, 256], 'text': '# of levels'},
            {'var': 'sym', 'value': True, 'type': 'bool', 'text': 'symmetric'},
            {'var': 'dist', 'value': 1, 'type': 'int', 'range': [1, 99], 'text': 'distance'},
            {'var': 'ang', 'value': 'all', 'type': 'list', 'range': ['0', '45', '90', '135', 'all'], 'text': 'direction [deg]'},
            {'var': 'property', 'value': 'contrast', 'type': 'list', 'range': ['contrast', 'dissimilarity', 'homogenity', 'ASM', 'energy', 'correlation'], 'text': 'texture property'}
            ]

        def __init__(self, *args, **kwargs):
            pyrat.FilterWorker.__init__(self, *args, **kwargs)
            self.blockprocess = True
            self.blockoverlap = self.win//2+1

        def pre(self, *args, **kwargs):
            amp = pyrat.data.getData()
            amp = np.uint8(np.float32(sarscale(amp)) / (256.0 / self.levels))
            newlayer = pyrat.data.addLayer(amp)
            pyrat.data.activateLayer(newlayer)

        def filter(self, array, *args, **kwargs):
            shp = array.shape
            array = array.reshape((np.prod(shp[0:-2]),) + shp[-2:])   # reshape to (nch, ny, nx)

            for k, amp in enumerate(np.abs(array)):  # loop over channels
                if self.ang == 'all':
                    angles = [0, np.pi/4, np.pi/2, 3*np.pi/4]
                else:
                    angles = [float(self.ang)/180.0*np.pi]

                shape = (amp.shape[0] - self.win + 1, amp.shape[1] - self.win + 1) + (self.win, self.win)
                strid = (amp.strides[0], amp.strides[1]) + amp.strides
                amp = as_strided(amp, shape=shape, strides=strid)
                for y in range(amp.shape[0]):
                    for x in range(amp.shape[1]):
                        win = amp[y, x]
                        p = greycomatrix(win, [self.dist], angles, levels=self.levels, symmetric=self.sym)  # calc glcm histogram
                        array[k, y+self.win//2, x+self.win//2] = np.mean(greycoprops(p, self.property))  # calc property
            array = array.reshape(shp)      # back to original shape
            return array

    def glcm(*args, **kwargs):
        return GLCM(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("GLCM module requires skimage library (install?)")

############################################################################################
############################################################################################
############################################################################################
