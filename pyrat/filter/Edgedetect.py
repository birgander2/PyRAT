from __future__ import print_function
import pyrat
import numpy as np
from scipy import ndimage


class Sobel(pyrat.FilterWorker):
    """
    Sobel edge detector

    Takes the square root of the sum of the squares of the horizontal and vertical Sobels
    to get a magnitude that is somewhat insensitive to direction.
    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Edge detection', 'entry': 'Sobel'}

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.blockprocess = True
        self.blockoverlap = 3

    def filter(self, array, *args, **kwargs):
        shp = array.shape
        array = array.reshape((np.prod(shp[0:-2]),) + shp[-2:])   # reshape to (nch, ny, nx)

        for k, arr in enumerate(np.abs(array)): # loop over channels
            dx = ndimage.sobel(arr, 0)  # horizontal derivative
            dy = ndimage.sobel(arr, 1)  # vertical derivative
            array[k, ...] = np.hypot(dx, dy)  # magnitude

        array = array.reshape(shp)      # back to original shape
        return array


def sobel(*args, **kwargs):
    return Sobel(*args, **kwargs).run(**kwargs)

############################################################################################
############################################################################################
############################################################################################

try:
    from skimage import filter

    class Roberts(pyrat.FilterWorker):
        """
        Roberts edge detector

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Edge detection', 'entry': 'Roberts'}

        def __init__(self, *args, **kwargs):
            pyrat.FilterWorker.__init__(self, *args, **kwargs)
            self.blockprocess = True
            self.blockoverlap = 3

        def filter(self, array, *args, **kwargs):
            shp = array.shape
            array = array.reshape((np.prod(shp[0:-2]),) + shp[-2:]) # reshape to (nch, ny, nx)

            for k, arr in enumerate(np.abs(array)):   # loop over channels
                array[k, ...] = filter.roberts(arr)

            array = array.reshape(shp)  # back to original shape
            return array

    def roberts(*args, **kwargs):
        return Roberts(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("Roberts module requires skimage library (install?)")


############################################################################################
############################################################################################
############################################################################################


class RoA(pyrat.FilterWorker):
    """
    Ratio of averages edge detector

    This filter uses 4 ratios: horizontal, vertical and the 2 diagonals

    For more information, see  A. Bovik: On detecting edges in speckled imagery,
    IEEE Transactions on Acoustics, Speech and Signal Processing, 36(10), pp. 1618-1627, 1988

    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Edge detection', 'entry': 'RoA'}
    para = [{'var': 'win', 'value': 3, 'type': 'int', 'range': [3, 999], 'text': 'Window size'}]

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.blockprocess = True
        self.blockoverlap = self.win

    def filter(self, array, *args, **kwargs):
        shp = array.shape
        array = array.reshape((np.prod(shp[0:-2]),) + shp[-2:])   # reshape to (nch, ny, nx)

        for k, arr in enumerate(np.abs(array)):  # loop over channels
            array[k, ...] = self.roa(arr, self.win)

        array = array.reshape(shp)      # back to original shape
        return array

    @staticmethod
    def roa(amp, win):
        shp = amp.shape
        dsh = (win+1)//2
        samp = ndimage.filters.uniform_filter(amp, size=win)
        grad = np.empty(((4, )+shp), dtype='f4')

        samp_p0 = np.roll(samp, dsh, axis=0)
        samp_0p = np.roll(samp, dsh, axis=1)
        samp_pp = np.roll(samp_0p, dsh, axis=0)
        samp_m0 = np.roll(samp, -dsh, axis=0)
        samp_0m = np.roll(samp, -dsh, axis=1)
        samp_mm = np.roll(samp_0m, -dsh, axis=0)
        samp_mp = np.roll(samp_0p, -dsh, axis=0)
        samp_pm = np.roll(samp_0m, dsh, axis=0)

        np.seterr(divide='ignore', invalid='ignore')
        grad[0, ...] = (samp_0p + samp_pp + samp_mp) / (samp_0m + samp_pm + samp_mm)  # vertical
        grad[1, ...] = (samp_pp + samp_p0 + samp_0p) / (samp_mm + samp_m0 + samp_0m)  # diagonal 1
        grad[2, ...] = (samp_p0 + samp_pp + samp_pm) / (samp_m0 + samp_mm + samp_mp)  # horizotal
        grad[3, ...] = (samp_mp + samp_0p + samp_m0) / (samp_pm + samp_p0 + samp_0m)  # diagonal 2
        grad = np.maximum(grad, 1.0 / grad)
        grad = np.sqrt(grad[0, ...]**2 + grad[1, ...]**2 + grad[2, ...]**2 + grad[3, ...]**2)
        np.seterr(divide='warn', invalid='warn')
        grad[~np.isfinite(grad)] = 0.0
        return grad


def roa(*args, **kwargs):
    return RoA(*args, **kwargs).run(**kwargs)

############################################################################################
############################################################################################
############################################################################################


class LeeRoA(pyrat.FilterWorker):
    """
    Ratio of averages edge detector in the version of J.S. Lee

    This filter uses 4 ratios: horizontal, vertical and the 2 diagonals

    For more information, see J.S. Lee et al.: Polarimetric SAR Speckle Filtering and Its
    Implication for Classification, IEEE Trans. Geos. Rem. Sens. 37(5), pp. 2363-2373, 1999

    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Edge detection', 'entry': 'Lee-RoA'}
    para = [{'var': 'win', 'value': 3, 'type': 'int', 'range': [3, 999], 'text': 'Window size'}]

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.blockprocess = True
        self.blockoverlap = self.win

    def filter(self, array, *args, **kwargs):
        shp = array.shape
        array = array.reshape((np.prod(shp[0:-2]),) + shp[-2:])   # reshape to (nch, ny, nx)

        for k, arr in enumerate(np.abs(array)):  # loop over channels
            array[k, ...] = self.leeroa(arr, self.win)

        array = array.reshape(shp)      # back to original shape
        return array

    @staticmethod
    def leeroa(amp, win):
        shp = amp.shape
        dsh = (win+1)//2
        samp = ndimage.filters.uniform_filter(amp, size=win)
        grad = np.empty(((8, )+shp), dtype='f4')
        xs = [dsh, dsh, 0, -dsh, -dsh, dsh, 0, -dsh]
        ys = [0, dsh, dsh, dsh, 0, -dsh, -dsh, -dsh]

        np.seterr(divide='ignore', invalid='ignore')
        for k in range(8):
            grad[k, ...] = np.abs(np.roll(np.roll(samp, ys[k], axis=0),xs[k], axis=1)/samp-1.0)
        grad = np.amax(np.abs(grad), axis=0)
        np.seterr(divide='warn', invalid='warn')
        grad[~np.isfinite(grad)] = 0.0
        return grad


def leeroa(*args, **kwargs):
    return LeeRoA(*args, **kwargs).run(**kwargs)
