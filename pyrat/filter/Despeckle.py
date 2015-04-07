from __future__ import print_function
import pyrat
import scipy as sp
import numpy as np
import logging
from scipy.ndimage import filters


class Boxcar(pyrat.FilterWorker):
    """
    Simple Boxcar filter...

    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Speckle filter', 'entry': 'Boxcar'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']},
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
            return filters.uniform_filter(array.real, win) + 1j * filters.uniform_filter(array.imag, win)
        elif self.phase is True:
            tmp = np.exp(1j * array)
            tmp = filters.uniform_filter(tmp.real, win) + 1j * filters.uniform_filter(tmp.imag, win)
            return np.angle(tmp)
        else:
            return filters.uniform_filter(array.real, win)


def boxcar(*args, **kwargs):
    return Boxcar(*args, **kwargs).run(**kwargs)

######################################################################################################################
######################################################################################################################
######################################################################################################################

class Lee(pyrat.FilterWorker):
    """
    Lee's classical speckle filter from 1981. Not the best one...

    :author: Andreas Reigber
    :param array: The image to filter (2D np.ndarray)
    :type array: float
    :param win: The filter window size (default: [7,7])
    :type win: integer
    :param looks=1.0: The effective number of looks of the input image.
    :type looks: float
    :returns: filtered image
    """
    gui = {'menu': 'SAR|Speckle filter', 'entry': 'Lee MMSE'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']},
        {'var': 'looks', 'value': 1.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'}
    ]

    def __init__(self, *args, **kwargs):
        super(Lee, self).__init__(*args, **kwargs)
        self.name = "LEE FILTER"
        self.blockprocess = True
        self.blockoverlap = self.win[0] // 2 + 1

    def filter(self, array, *args, **kwargs):
        array[np.isnan(array)] = 0.0
        shape = array.shape
        if len(shape) == 3:
            array = np.abs(array)
            span = np.sum(array ** 2, axis=0)
            array = array[np.newaxis, ...]
        elif len(shape) == 4:
            span = np.trace(array, axis1=0, axis2=1)
        else:
            array = np.abs(array)
            span = array ** 2
            array = array[np.newaxis, np.newaxis, ...]
        lshape = array.shape[0:2]

        sig2 = 1.0 / self.looks
        sfak = 1.0 + sig2
        m2arr = sp.ndimage.filters.uniform_filter(span ** 2, size=self.win)
        marr = sp.ndimage.filters.uniform_filter(span, size=self.win)
        vary = (m2arr - marr ** 2).clip(1e-10)
        varx = ((vary - marr ** 2 * sig2) / sfak).clip(0)
        kfac = varx / vary

        out = np.empty_like(array)
        for k in range(lshape[0]):
            for l in range(lshape[1]):
                if np.iscomplexobj(array):
                    out[k, l, ...] = sp.ndimage.filters.uniform_filter(array[k, l, ...].real, size=self.win) + \
                                     1j * sp.ndimage.filters.uniform_filter(array[k, l, ...].imag, size=self.win)
                else:
                    out[k, l, ...] = sp.ndimage.filters.uniform_filter(array[k, l, ...], size=self.win)
                out[k, l, ...] += (array[k, l, ...] - out[k, l, ...]) * kfac
        return np.squeeze(out)


def lee(*args, **kwargs):
    return Lee(*args, **kwargs).run(**kwargs)


class Gauss(pyrat.FilterWorker):
    """
    Simple Gaussian filter...

    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Speckle filter', 'entry': 'Gauss'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Sigma', 'subtext': ['range', 'azimuth']},
        {'var': 'phase', 'value': False, 'type': 'bool', 'text': 'Phase'}
        ]

    def __init__(self, *args, **kwargs):
        super(Gauss, self).__init__(*args, **kwargs)
        self.name = "GAUSS FILTER"
        self.blockprocess = True
        self.blockoverlap = 2*self.win[0] + 1

    def filter(self, array, *args, **kwargs):
        win = self.win
        if array.ndim == 3:
            win = [1] + self.win
        if array.ndim == 4:
            win = [1, 1] + self.win
        array[np.isnan(array)] = 0.0
        if np.iscomplexobj(array):
            return filters.gaussian_filter(array.real, win) + 1j * filters.gaussian_filter(array.imag, win)
        elif self.phase is True:
            tmp = np.exp(1j * array)
            tmp = filters.gaussian_filter(tmp.real, win) + 1j * filters.gaussian_filter(tmp.imag, win)
            return np.angle(tmp)
        else:
            return filters.gaussian_filter(array.real, win)


def gauss(*args, **kwargs):
    return Gauss(*args, **kwargs).run(**kwargs)


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


class RefinedLee(pyrat.FilterWorker):
    """
    Refined Lee speckle filter

    further information:
    J.S.Lee et al.: "Speckle Reduction in Multipolarization Multifrequency SAR Imagery",
    Trans. on Geoscience and Remote Sensing, Vol. 29, No. 4, pp. 535-544, 1991

    :author: Andreas Reigber
    :param array: The image to filter (2D np.ndarray)
    :type array: float
    :param win: The filter window size (default: [7,7])
    :type win: integer
    :param looks=1.0: The effective number of looks of the input image.
    :type looks: float
    :param threshold=0.5: Threshold on which switch to normal Lee filtering
    :type threshold: float
    :param method='original': The edge detection method used.
    :type method: string
    :returns: filtered image
    """
    gui = {'menu': 'SAR|Speckle filter', 'entry': 'Refined Lee'}
    para = [
        {'var': 'win', 'value': 7, 'type': 'int', 'range': [3, 999], 'text': 'Window size'},
        {'var': 'looks', 'value': 1.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
        {'var': 'threshold', 'value': 0.5, 'type': 'float', 'range': [0.0, 9.0], 'text': 'threshold'},
        {'var': 'method', 'value': 'original', 'type': 'list', 'range': ['original', 'cov'], 'text': 'edge detector'}
    ]

    def __init__(self, *args, **kwargs):
        super(RefinedLee, self).__init__(*args, **kwargs)
        self.name = "REFINED LEE FILTER"
        self.blockprocess = True
        self.blockoverlap = self.win // 2 + 1

    def filter(self, array, *args, **kwargs):
        array[np.isnan(array)] = 0.0
        shape = array.shape
        if len(shape) == 3:
            array = np.abs(array)
            span = np.sum(array ** 2, axis=0)
            array = array[np.newaxis, ...]
        elif len(shape) == 4:
            span = np.trace(array, axis1=0, axis2=1)
        else:
            array = np.abs(array)
            span = array ** 2
            array = array[np.newaxis, np.newaxis, ...]
        lshape = array.shape[0:2]

        # ---------------------------------------------
        # INIT & SPAN
        # ---------------------------------------------

        sig2 = 1.0 / self.looks
        sfak = 1.0 + sig2

        # nrx = array.shape[-1]
        #
        # lshape = array.shape[0:-2]
        # if len(lshape) == 2:
        # # span = np.abs(np.trace(array,axis1=0,axis2=1))
        #     span = np.abs(array[0, 0, ...] + array[1, 1, ...] + array[2, 2, ...])
        # else:
        #     logging.error("Data not in matrix form")

        # ---------------------------------------------
        # TURNING BOX
        # ---------------------------------------------

        cbox = np.zeros((9, self.win, self.win), dtype='float32')
        chbox = np.zeros((self.win, self.win), dtype='float32')
        chbox[0:self.win // 2 + 1, :] = 1
        cvbox = np.zeros((self.win, self.win), dtype='float32')
        for k in range(self.win):
            cvbox[k, 0:k + 1] = 1

        cbox[0, ...] = np.rot90(chbox, 3)
        cbox[1, ...] = np.rot90(cvbox, 1)
        cbox[2, ...] = np.rot90(chbox, 2)
        cbox[3, ...] = np.rot90(cvbox, 0)
        cbox[4, ...] = np.rot90(chbox, 1)
        cbox[5, ...] = np.rot90(cvbox, 3)
        cbox[6, ...] = np.rot90(chbox, 0)
        cbox[7, ...] = np.rot90(cvbox, 2)
        for k in range(self.win // 2 + 1):
            for l in range(self.win // 2 - k, self.win // 2 + k + 1):
                cbox[8, k:self.win - k, l] = 1

        for k in range(9):
            cbox[k, ...] /= np.sum(cbox[k, ...])

        ampf1 = np.empty((9,) + span.shape)
        ampf2 = np.empty((9,) + span.shape)
        for k in range(9):
            ampf1[k, ...] = filters.correlate(span ** 2, cbox[k, ...])
            ampf2[k, ...] = filters.correlate(span, cbox[k, ...]) ** 2

        # ---------------------------------------------
        # GRADIENT ESTIMATION
        # ---------------------------------------------
        np.seterr(divide='ignore', invalid='ignore')

        if self.method == 'original':
            xs = [+2, +2, 0, -2, -2, -2, 0, +2]
            ys = [0, +2, +2, +2, 0, -2, -2, -2]
            samp = filters.uniform_filter(span, self.win // 2)
            grad = np.empty((8,) + span.shape)
            for k in range(8):
                grad[k, ...] = np.abs(np.roll(np.roll(samp, ys[k], axis=0), xs[k], axis=1) / samp - 1.0)
            magni = np.amax(grad, axis=0)
            direc = np.argmax(grad, axis=0)
            direc[magni < self.threshold] = 8
        elif self.method == 'cov':
            grad = np.empty((8,) + span.shape)
            for k in range(8):
                grad[k, ...] = np.abs((ampf1[k, ...] - ampf2[k, ...]) / ampf2[k, ...])
                direc = np.argmin(grad, axis=0)
        else:
            logging.error("Unknown method!")

        np.seterr(divide='warn', invalid='warn')
        # ---------------------------------------------
        # FILTERING
        # ---------------------------------------------
        out = np.empty_like(array)
        dbox = np.zeros((1, 1) + (self.win, self.win))
        for l in range(9):
            grad = ampf1[l, ...]
            mamp = ampf2[l, ...]
            dbox[0, 0, ...] = cbox[l, ...]

            vary = (grad - mamp).clip(1e-10)
            varx = ((vary - mamp * sig2) / sfak).clip(0)
            kfac = varx / vary
            if np.iscomplexobj(array):
                mamp = filters.correlate(array.real, dbox) + 1j * filters.convolve(array.imag, dbox)
            else:
                mamp = filters.correlate(array, dbox)
            idx = np.argwhere(direc == l)
            out[:, :, idx[:, 0], idx[:, 1]] = (mamp + (array - mamp) * kfac)[:, :, idx[:, 0], idx[:, 1]]

        return out


def reflee(*args, **kwargs):
    return RefinedLee(*args, **kwargs).run(**kwargs)


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
try:
    from .Cy_Despeckle import cy_leesigma
    class LeeSigma(pyrat.FilterWorker):
        """
        Lee's original sigma speckle filter...
        J.S. Lee: A Simple Speckle Smoothing  Algorithm for Synthetic Aperture Radar Images.
        IEEE Transactions on System, Man, and Cybernetics, Vol. SMC-13, No.1, pp.85-89, 1983']

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'Lee Sigma (old)'}
        para = [
            {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']},
            {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude','intensity'], 'text': 'SAR data type'}
            ]

        def __init__(self, *args, **kwargs):
            super(LeeSigma, self).__init__(*args, **kwargs)
            self.name = "SIGMA LEE SPECKLE FILTER"
            self.blockprocess = True
            self.blockoverlap = self.win[0] // 2 + 1

        def filter(self, array, *args, **kwargs):
            np.seterr(divide='ignore', invalid='ignore')
            if np.iscomplexobj(array):
                array = np.abs(array)
                self.type = "amplitude"
            if self.type == "amplitude":
                array *= array
            array = cy_leesigma(array, looks=self.looks, win=self.win)
            if self.type == "amplitude":
                array = np.sqrt(array)
            array[~np.isfinite(array)] = 0.0
            np.seterr(divide='warn', invalid='warn')

            return array


    def leesigma(*args, **kwargs):
        return LeeSigma(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("LeeSigma module not found. (run build process?)")




#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# from .Cy_Despeckle import cy_leesigmanew
# class LeeSigmaNew(pyrat.FilterWorker):
#     """
#     Lee's original sigma speckle filter...
#     J.S. Lee: A Simple Speckle Smoothing  Algorithm for Synthetic Aperture Radar Images.
#     IEEE Transactions on System, Man, and Cybernetics, Vol. SMC-13, No.1, pp.85-89, 1983']
#
#     :author: Andreas Reigber
#     """
#     gui = {'menu': 'SAR|Speckle filter', 'entry': 'Lee Sigma (new)'}
#     para = { 'win': {'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
#                      'subtext': ['range', 'azimuth']},
#              'looks': {'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
#              'type': {'value': 'amplitude', 'type': 'list', 'range': ['amplitude','intensity'], 'text': 'SAR data type'}}
#
#     def __init__(self, *args, **kwargs):
#         super(LeeSigmaNew, self).__init__(*args, **kwargs)
#         self.name = "SIGMA LEE SPECKLE FILTER"
#         self.blockprocess = True
#         self.blockoverlap = self.win[0] // 2 + 1
#
#     def filter(self, array, *args, **kwargs):
#         oarray = array.copy()
#         if self.type == "amplitude":
#             array *= array
#         array = cy_leesigmanew(array, looks=self.looks, win=self.win)
#         if self.type == "amplitude":
#             array = np.sqrt(array)
#         array[~np.isfinite(array)] = 0.0
#
#         return array
#
#
# def leesigmanew(*args, **kwargs):
#     return LeeSigmaNew(*args, **kwargs).run(**kwargs)


