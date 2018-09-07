import logging

import pyrat
import scipy as sp
from scipy import optimize as opt
import numpy as np


class Boxcar(pyrat.FilterWorker):
    """
    Boxcar / Moving average (speckle) filter.

    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Speckle filter', 'entry': 'Boxcar'}
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
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
         'subtext': ['range', 'azimuth']},
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
            span = np.abs(np.trace(array, axis1=0, axis2=1))
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


@pyrat.docstringfrom(Lee)
def lee(*args, **kwargs):
    return Lee(*args, **kwargs).run(**kwargs)


class Kuan(pyrat.FilterWorker):
    """
    Kuan filter. Similar performance to the Lee filter but with a different weighting function.

    further information:
    D.T.Kuan et al. “Adaptive restoration of images with speckle,” JEEE Trans. Acoustics, Speech and Sig.
    Proc., vol. ASSP-35, pp. 373-383, March 1987

    :author: Joel Amao
    :param array: The image to filter (2D np.ndarray)
    :type array: float
    :param win: The filter window size (default: [7,7])
    :type win: integer
    :param looks=1.0: The effective number of looks of the input image.
    :type looks: float
    :returns: filtered image
    """
    gui = {'menu': 'SAR|Speckle filter', 'entry': 'Kuan'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
         'subtext': ['range', 'azimuth']},
        {'var': 'looks', 'value': 1.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'}
    ]

    def __init__(self, *args, **kwargs):
        super(Kuan, self).__init__(*args, **kwargs)
        self.name = "KUAN FILTER"
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
            span = np.abs(np.trace(array, axis1=0, axis2=1))
        else:
            array = np.abs(array)
            span = array ** 2
            array = array[np.newaxis, np.newaxis, ...]
        lshape = array.shape[0:2]

        # Calculates the squared mean of the array (marr) and the mean squared array (m2arr)
        m2arr = sp.ndimage.filters.uniform_filter(span ** 2, size=self.win)
        marr = sp.ndimage.filters.uniform_filter(span, size=self.win)
        # Variance within window, follows the identity Var(x) = E(x**2) - [E(X)]**2
        vary = (m2arr - marr ** 2).clip(1e-10)
        # Standard deviation within window
        stdDev = np.sqrt(vary)
        # cu and ci are the main parameters of the weight function w
        cu = np.sqrt(1 / self.looks)
        if not cu:
            cu = 0.25

        ci = (stdDev / marr)
        if not ci.any():
            ci = 0.0001

        # Clipped weighted function
        w = (1 - ((cu ** 2) / (ci ** 2))).clip(0) / ((1 + (cu ** 2)).clip(1e-10))

        # Filters each channel separately
        out = np.empty_like(array)
        for k in range(lshape[0]):
            for l in range(lshape[1]):
                if np.iscomplexobj(array):
                    out[k, l, ...] = sp.ndimage.filters.uniform_filter(array[k, l, ...].real, size=self.win) + \
                                     1j * sp.ndimage.filters.uniform_filter(array[k, l, ...].imag, size=self.win)
                else:
                    out[k, l, ...] = sp.ndimage.filters.uniform_filter(array[k, l, ...], size=self.win)
                # Main output
                out[k, l, ...] = (array[k, l, ...] * w) + out[k, l, ...] * (1 - w)
        return np.squeeze(out)


@pyrat.docstringfrom(Kuan)
def kuan(*args, **kwargs):
    return Kuan(*args, **kwargs).run(*args, **kwargs)


class IDAN(pyrat.FilterWorker):
    """
    IDAN filter. Intensity driven adaptive neighborhood filter.

    further information:
    G. Vasile et al. “Intensity-Driven Adaptive-Neighborhood Technique for Polarimetric and Interferometric
    SAR Parameters Estimation” IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING,Vol 44. No. 6, 2006


    :author: Joel Amao
    :param array: The image to filter (2D np.ndarray)
    :type array: float
    :param looks = 1.0: The effective number of looks of the input image.
    :param nmax = 50: The maximum number of neighbours to add
    :param LLMMSE: Activates LLMMSE processing after the region growing
    :type looks: float
    :returns: filtered image
    """
    gui = {'menu': 'SAR|Speckle filter', 'entry': 'IDAN filter'}
    para = [
        {'var': 'looks', 'value': 1.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
        {'var': 'nmax', 'value': 50, 'type': 'int', 'range': [40, 200], 'text': 'Max # of Neighbors'},
        {'var': 'LLMMSE', 'value': True, 'type': 'bool', 'text': 'Activates LLMMSE processing'}
    ]

    def __init__(self, *args, **kwargs):
        super(IDAN, self).__init__(*args, **kwargs)
        self.name = "IDAN FILTER"
        self.blockprocess = True
        self.blocksize = self.nmax * 4
        self.blockoverlap = self.nmax
        # self.nthreads = 1
        # todo: add a mininum nmax parameters
        # todo: test InSAR data

    def filter(self, array, *args, **kwargs):

        MMSEproc = self.LLMMSE
        array[np.isnan(array)] = 0.0
        shape = array.shape
        if len(shape) == 3:
            array = np.abs(array)
            span = np.sum(array ** 2, axis=0)
            array = array[np.newaxis, ...]
        elif len(shape) == 4:
            span = np.abs(np.trace(array, axis1=0, axis2=1))
        else:
            array = np.abs(array)
            span = array ** 2
            array = array[np.newaxis, np.newaxis, ...]

        win = np.array([3, 3])
        ldim = array.shape[0:2]
        rdim = array.shape[2:4]
        nmax = self.nmax
        nlook = 1 / (np.sqrt(self.looks)) / 3.0

        # ==============================================================================
        # 1.1 Rough estimation of the seed value
        # ==============================================================================

        intensities = np.zeros((ldim[0], rdim[0], rdim[1]), dtype=np.float32)
        med_arr = np.zeros_like(intensities)
        for i in range(ldim[0]):
            intensities[i, :] = array[i, i, :, :].real
            med_arr[i, :] = sp.ndimage.filters.median_filter(intensities[i, :], size=win)

        data = array
        ldim = data.shape[0:2]
        rdim = data.shape[2:4]
        maxp = rdim[0] * rdim[1]

        mask_vec = np.zeros(shape=[maxp, nmax], dtype=np.int64)
        background = np.zeros_like(mask_vec)
        nnList = np.zeros((maxp, 1), dtype=np.int64)

        # k = (xx)*rdim[1] + (yy)

        # ==============================================================================
        # 1.2 Region growing
        # ==============================================================================

        thres = 2 * nlook
        intensities = intensities.reshape((ldim[0], maxp))
        med_arr = med_arr.reshape(ldim[0], maxp)
        neighbours = [1, -1, -rdim[1], rdim[1], rdim[1] + 1, rdim[1] - 1, -rdim[1] + 1, -rdim[1] - 1]

        for k in range(maxp):
            stack = [(k)]
            nn = 0
            bg_nn = 0
            while stack:
                y = stack.pop(0)
                nnList[k] = nn
                if nn == 0:
                    mask_vec[k, 0] = y

                for dy in neighbours:
                    ny = y + dy
                    # and >>> and
                    if (not (ny in background[k, :]) and 0 <= ny and not (
                        ny in mask_vec[k, :]) and nn < nmax and ny < maxp):

                        up = np.abs(intensities[:, ny] - med_arr[:, k])
                        down = np.maximum(np.abs(med_arr[:, k]), 1e-10)
                        sum = (up / down).sum()

                        if (sum <= thres):
                            stack.append(ny)
                            if nn < nmax:
                                nn += 1
                            if nn < nmax:
                                mask_vec[k, nn] = ny
                        else:
                            if bg_nn < nmax:
                                background[k, bg_nn] = ny
                                bg_nn += 1

        # ==============================================================================
        # 2.1 Refined estimation of the seed value
        # ==============================================================================
        dataV = data.reshape(ldim[0], ldim[1], maxp)
        updatedSeed = np.zeros_like(med_arr)

        for k in range(maxp):
            for i in range(ldim[0]):
                if np.where(mask_vec[k, :] > 0)[0].size > 0:
                    updatedSeed[i, k] = np.mean(dataV[i, i, mask_vec[k, :][np.where(mask_vec[k, :] > 0)]].real)

        # ==============================================================================
        # 2.2 Reinspection of the background pixels
        # ==============================================================================
        thres = 6 * nlook

        for k in range(maxp):
            if nnList[k] < nmax:
                stack = background[k][np.where(background[k] > 0)].tolist()
                while stack:
                    nnList[k] += 1
                    y = stack.pop(0)
                    up = np.abs(intensities[:, y] - updatedSeed[:, k])
                    down = np.maximum(np.abs(updatedSeed[:, k]), 1e-10)
                    sum = (up / down).sum()
                    if (sum <= thres):
                        if nnList[k] < nmax:
                            mask_vec[k, nnList[k]] = y

        # ==============================================================================
        # 3.1 Parameter estimation
        # ==============================================================================
        avg = np.zeros(shape=[ldim[0], ldim[1], maxp], dtype=np.complex64)

        for k in range(maxp):
            for i in range(ldim[0]):
                for j in range(ldim[1]):
                    if np.where(mask_vec[k, :] > 0)[0].size > 0:
                        idx = mask_vec[k, :][np.where(mask_vec[k, :] > 0)]
                        avg[i, j, k] = np.mean(dataV[i, j, idx])
                    else:
                        avg[i, j, k] = dataV[i, j, k]

        # Calculates the squared mean of the array (marr) and the mean squared array (m2arr)
        m2arr = sp.ndimage.filters.uniform_filter(span ** 2, size=win)
        marr = sp.ndimage.filters.uniform_filter(span, size=win)
        # Variance within window, follows the identity Var(x) = E(x**2) - [E(X)]**2
        vary = (m2arr - marr ** 2).clip(1e-10)
        # Standard deviation within window
        stdDev = np.sqrt(vary)
        # cu and ci are the main parameters of the weight function w
        cu = np.sqrt(1 / self.looks)
        ci = (stdDev / marr)

        # Clipped weighted function
        w = (1 - ((cu ** 2) / (ci ** 2))).clip(0) / ((1 + (cu ** 2)).clip(1e-10))
        w = w.reshape(maxp)

        # LLMMSE
        if MMSEproc:
            LLMMSE = np.zeros_like(avg, dtype=np.complex64)
            for k in range(maxp):
                for i in range(ldim[0]):
                    for j in range(ldim[1]):
                        LLMMSE[i, j, k] = avg[i, j, k] + w[k] * (dataV[i, j, k] - avg[i, j, k])
            LLMMSE = LLMMSE.reshape(ldim[0], ldim[1], rdim[0], rdim[1])
            return np.squeeze(LLMMSE)
        else:
            # Reshaping the array and removing the borders
            avg = avg.reshape(ldim[0], ldim[1], rdim[0], rdim[1])
            return np.squeeze(avg)


@pyrat.docstringfrom(IDAN)
def idan(*args, **kwargs):
    return IDAN(*args, **kwargs).run(*args, **kwargs)


class Gauss(pyrat.FilterWorker):
    """
    Gaussian (speckle) filter...

    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Speckle filter', 'entry': 'Gauss'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Sigma',
         'subtext': ['range', 'azimuth']},
        {'var': 'phase', 'value': False, 'type': 'bool', 'text': 'Phase'}
    ]

    def __init__(self, *args, **kwargs):
        super(Gauss, self).__init__(*args, **kwargs)
        self.name = "GAUSS FILTER"
        self.blockprocess = True
        self.blockoverlap = 2 * self.win[0] + 1

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
            span = np.abs(np.trace(array, axis1=0, axis2=1))
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
            ampf1[k, ...] = sp.ndimage.filters.correlate(span ** 2, cbox[k, ...])
            ampf2[k, ...] = sp.ndimage.filters.correlate(span, cbox[k, ...]) ** 2

        # ---------------------------------------------
        # GRADIENT ESTIMATION
        # ---------------------------------------------
        np.seterr(divide='ignore', invalid='ignore')

        if self.method == 'original':
            xs = [+2, +2, 0, -2, -2, -2, 0, +2]
            ys = [0, +2, +2, +2, 0, -2, -2, -2]
            samp = sp.ndimage.filters.uniform_filter(span, self.win // 2)
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
                mamp = sp.ndimage.filters.correlate(array.real, dbox) + 1j * sp.ndimage.filters.convolve(array.imag,
                                                                                                         dbox)
            else:
                mamp = sp.ndimage.filters.correlate(array, dbox)
            idx = np.argwhere(direc == l)
            out[:, :, idx[:, 0], idx[:, 1]] = (mamp + (array - mamp) * kfac)[:, :, idx[:, 0], idx[:, 1]]

        return out


@pyrat.docstringfrom(RefinedLee)
def refinedlee(*args, **kwargs):
    return RefinedLee(*args, **kwargs).run(**kwargs)


# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------
# try:
#     from .Despeckle_extensions import cy_leesigmaold
#     class LeeSigmaOld(pyrat.FilterWorker):
#         """
#         Lee's original sigma speckle filter...
#         J.S. Lee: A Simple Speckle Smoothing  Algorithm for Synthetic Aperture Radar Images.
#         IEEE Transactions on System, Man, and Cybernetics, Vol. SMC-13, No.1, pp.85-89, 1983']
#
#         :author: Andreas Reigber
#         """
#         gui = {'menu': 'SAR|Speckle filter', 'entry': 'Lee Sigma (old)'}
#         para = [
#             {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']},
#             {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
#             {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude','intensity'], 'text': 'SAR data type'}
#             ]
#
#         def __init__(self, *args, **kwargs):
#             super(LeeSigmaOld, self).__init__(*args, **kwargs)
#             self.name = "SIGMA LEE SPECKLE FILTER"
#             self.blockprocess = True
#             self.blockoverlap = self.win[0] // 2 + 1
#
#         def filter(self, array, *args, **kwargs):
#             np.seterr(divide='ignore', invalid='ignore')
#             if np.iscomplexobj(array):
#                 array = np.abs(array)
#                 self.type = "amplitude"
#             if self.type == "amplitude":
#                 array *= array
#             array = cy_leesigmaold(array, looks=self.looks, win=self.win)
#             if self.type == "amplitude":
#                 array = np.sqrt(array)
#             array[~np.isfinite(array)] = 0.0
#             np.seterr(divide='warn', invalid='warn')
#
#             return array
#
#
#     def leesigmaold(*args, **kwargs):
#         return LeeSigmaOld(*args, **kwargs).run(**kwargs)
#
# except ImportError:
#     logging.info("LeeSigma module not found. (run build process?)")


# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------

try:
    from .Despeckle_extensions import cy_leesigma


    class LeeSigma(pyrat.FilterWorker):
        """
        Lee's original sigma speckle filter. Fast implementation in Cython.

        J.S. Lee: A Simple Speckle Smoothing  Algorithm for Synthetic Aperture Radar Images.
        IEEE Transactions on System, Man, and Cybernetics, Vol. SMC-13, No.1, pp.85-89, 1983']

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'Lee Sigma (old)'}
        para = [
            {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
             'subtext': ['range', 'azimuth']},
            {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity'],
             'text': 'SAR data type'}
        ]

        def __init__(self, *args, **kwargs):
            super(LeeSigma, self).__init__(*args, **kwargs)
            self.name = "SIGMA LEE SPECKLE FILTER"
            self.blockprocess = True
            self.blockoverlap = self.win[0] // 2 + 1

        def filter(self, array, *args, **kwargs):
            if array.ndim == 3:  # polarimetric vector
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array) ** 2
                    self.type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                self.type = "amplitude"
            elif array.ndim == 4:  # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                self.type = "intensity"
            else:  # single channel data
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array) ** 2
                    self.type = "amplitude"
                else:
                    array = np.abs(array)
                span = array.copy()
                array = array[np.newaxis, np.newaxis, ...]
            array = cy_leesigma(span, array, looks=self.looks, win=self.win)
            array[~np.isfinite(array)] = 0.0
            if self.type == "amplitude":
                array[array < 0] = 0.0
                array = np.sqrt(array)
            return np.squeeze(array)


    @pyrat.docstringfrom(LeeSigma)
    def leesigma(*args, **kwargs):
        return LeeSigma(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("LeeSigma cython module not found. (run build process?)")

# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------

try:
    from .Despeckle_extensions import cy_leeimproved


    class LeeSigma2(pyrat.Worker):
        """
        Lee's improved sigma speckle filter. Fast implementation in Cython.

        J.S. Lee et al.: Improved Sigma Filter for Speckle Filtering of SAR Imagery.
        IEEE Transactions on Geoscience and Remote Sensing, Vol. 47, No.1, pp. 202-213, 2009']

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'Lee Sigma Improved'}
        para = [
            {'var': 'win', 'value': [9, 9], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
             'subtext': ['range', 'azimuth']},
            {'var': 'looks', 'value': 1.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'sigma', 'value': 0.9, 'type': 'float', 'range': [0.0, 1.0], 'text': 'sigma range (0.0-1.0)'},
            {'var': 'perc', 'value': 0.02, 'type': 'float', 'range': [0.0, 1.0],
             'text': 'point target percentile (0.0-1.0)'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity'],
             'text': 'SAR data type'}]

        def __init__(self, *args, **kwargs):
            super(LeeSigma2, self).__init__(*args, **kwargs)
            self.name = "IMPROVED SIGMA LEE"
            self.blockprocess = True
            self.blockoverlap = self.win[0] // 2 + 1
            # self.nthreads = 1

        def run(self, *args, **kwargs):
            P = ProgressBar('  ' + self.name, 10)
            bounds = opt.fmin(self.optf, [0.5, 2.0], args=(self.looks, self.sigma), disp=False)  # calc sigma bounds
            newsig = self.newsig(bounds[0], bounds[1], sigrng=self.sigma, looks=self.looks)  # calc new stddev
            P.update(0)
            perc = 100.0 - self.perc * 100.0  # point target theshold
            pthreshold = np.mean(self.layer_accumulate(self.estimate_percentile, type=self.type, perc=perc))
            P.update(2)

            layer = self.layer_process(self.leeimproved, bounds=bounds, newsig=newsig,
                                       thres=pthreshold, looks=self.looks, win=self.win, type=self.type)
            P.update(10)
            del P
            pyrat.activate(layer)
            return layer

        @staticmethod
        def leeimproved(array, bounds=(1, 2), thres=10.0, looks=1.0, win=(7, 7), newsig=0.5, type='amplitude',
                        **kwargs):

            if array.ndim == 3:  # polarimetric vector
                if np.iscomplexobj(array) or type == "amplitude":
                    array = np.abs(array) ** 2
                    type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                type = "amplitude"
            elif array.ndim == 4:  # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                type = "intensity"
            else:  # single channel data
                if np.iscomplexobj(array) or type == "amplitude":
                    array = np.abs(array) ** 2
                    type = "amplitude"
                else:
                    array = np.abs(array)
                span = array.copy()
                array = array[np.newaxis, np.newaxis, ...]

            array = cy_leeimproved(span, array, bounds=bounds, thres=thres, looks=looks, win=win, newsig=newsig)
            array[~np.isfinite(array)] = 0.0
            if type == "amplitude":
                array[array < 0] = 0.0
                array = np.sqrt(array)
            return np.squeeze(array)

        def estimate_percentile(self, array, perc=98.0, type='amplitude', **kwargs):
            if array.ndim == 3:  # polarimetric vector
                if np.iscomplexobj(array) or type == "amplitude":
                    span = np.sum(np.abs(array) ** 2, axis=0)
                else:
                    span = np.sum(np.abs(array), axis=0)
            elif array.ndim == 4:  # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
            else:  # single channel data
                if np.iscomplexobj(array) or type == "amplitude":
                    span = np.abs(array) ** 2
                else:
                    span = np.abs(array)
            return np.percentile(span, perc)

        @staticmethod
        def specklepdf(i, looks=1.0):
            if i < 0.0:
                return 0.0
            else:
                return ((looks ** looks) * (i ** (looks - 1.0))) / sp.special.gamma(looks) * np.exp(-looks * i)

        @staticmethod
        def meanpdf(i, looks=1.0):
            if i < 0.0:
                return 0.0
            else:
                return ((looks ** looks) * (i ** (looks - 1.0))) / sp.special.gamma(looks) * np.exp(-looks * i) * i

        def sigpdf(self, i, looks=1.0):
            return (i - 1.0) ** 2 * self.specklepdf(i, looks=1.0)

        def newsig(self, i1, i2, sigrng=0.9, looks=1.0):
            return 1 / sigrng * sp.integrate.quad(self.sigpdf, i1, i2, args=(looks,))[0]

        def sigmarange(self, i1, i2, looks=1.0):
            return np.clip(sp.integrate.quad(self.specklepdf, i1, i2, args=(looks,))[0], 1e-10, 1.0)

        def intmean(self, i1, i2, looks=1.0):
            return 1.0 / self.sigmarange(i1, i2, looks) * sp.integrate.quad(self.meanpdf, i1, i2, args=(looks,))[0]

        def optf(self, i, looks, sigr):
            return (self.sigmarange(i[0], i[1], looks) - sigr) ** 2 + (self.intmean(i[0], i[1], looks) - 1.0) ** 2


    @pyrat.docstringfrom(LeeSigma2)
    def leesigma2(*args, **kwargs):
        return LeeSigma2(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("LeeSigma2 cython module not found. (run build process?)")

# # ---------------------------------------------------------------------------------------------------------------
# # ---------------------------------------------------------------------------------------------------------------
# # ---------------------------------------------------------------------------------------------------------------
#
# try:
#     from .Despeckle_extensions import cy_leeimproved_old
#
#     class LeeSigmaImprovedOld(pyrat.Worker):
#         """
#         Lee's improved sigma speckle filter...
#         J.S. Lee et al.: Improved Sigma Filter for Speckle Filtering of SAR Imagery.
#         IEEE Transactions on Geoscience and Remote Sensing, Vol. 47, No.1, pp. 202-213, 2009']
#
#         :author: Andreas Reigber
#         """
#         gui = {'menu': 'SAR|Speckle filter', 'entry': 'Lee Sigma Improved'}
#         para = [
#             {'var': 'win', 'value': [9, 9], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']},
#             {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
#             {'var': 'sigma', 'value': 0.9, 'type': 'float', 'range': [0.0, 1.0], 'text': 'sigma range (0.0-1.0)'},
#             {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude','intensity'], 'text': 'SAR data type'}
#             ]
#
#         def __init__(self, *args, **kwargs):
#             super(LeeSigmaImprovedOld, self).__init__(*args, **kwargs)
#             self.name = "IMPROVED SIGMA LEE"
#             self.blockprocess = True
#             self.blockoverlap = self.win[0] // 2 + 1
#             # self.nthreads = 1
#
#         def run(self, *args, **kwargs):
#             bounds = opt.fmin(self.optf, [0.5, 2.0], args=(self.looks, self.sigma), disp=False)
#
#             newsig = self.newsig(bounds[0], bounds[1], sigrng=self.sigma, looks=self.looks)
#             pthreshold = np.mean(self.layer_accumulate(self.estimate_percentile))
#
#             print(bounds, pthreshold, self.sigma)
#
#             layer = self.layer_process(self.leeimproved, bounds=bounds, newsig=newsig, thres=pthreshold, looks=self.looks, win=self.win)
#             pyrat.activate(layer)
#             return layer
#
#         @staticmethod
#         def leeimproved(data, bounds=(1, 2), thres=10.0, looks=1.0, win=(7, 7), newsig=0.5, **kwargs):
#             data **= 2
#             out = cy_leeimproved_old(data, bounds=bounds, thres=thres, looks=looks, win=win, newsig=newsig)
#             out **= 0.5
#             data[~np.isfinite(data)] = 0.0
#             return out
#
#         @staticmethod
#         def estimate_percentile(data, perc=98, **kwargs):
#             data **= 2
#             return np.percentile(data, perc)
#
#         @staticmethod
#         def specklepdf(i, looks=1.0):
#             return ((looks**looks) * (i**(looks-1.0))) / np.math.factorial(looks-1.0) * np.exp(-looks*i)
#
#         @staticmethod
#         def meanpdf(i, looks=1.0):
#             return ((looks**looks) * (i**(looks-1.0))) / np.math.factorial(looks-1.0) * np.exp(-looks*i) * i
#
#         def sigpdf(self, i, looks=1.0):
#             return (i - 1.0)**2 * self.specklepdf(i, looks=1.0)
#
#         def newsig(self, i1, i2, sigrng=0.9, looks=1.0):
#             return 1 / sigrng * sp.integrate.quad(self.sigpdf, i1, i2, args=(looks, ))[0]
#
#         def sigmarange(self, i1, i2, looks=1.0):
#             return np.clip(sp.integrate.quad(self.specklepdf, i1, i2, args=(looks, ))[0], 1e-10, 1.0)
#
#         def intmean(self, i1, i2, looks=1.0):
#             return 1.0 / self.sigmarange(i1, i2, looks) * sp.integrate.quad(self.meanpdf, i1, i2, args=(looks, ))[0]
#
#         def optf(self, i, looks, sigr):
#             return (self.sigmarange(i[0], i[1], looks) - sigr)**2 + (self.intmean(i[0], i[1], looks) - 1.0)**2
#
#     def leesigmaimpX(*args, **kwargs):
#         return LeeSigmaImprovedOld(*args, **kwargs).run(**kwargs)
#
# except ImportError:
#     logging.info("LeeSigmaImproved cython module not found. (run build process?)")
#
#
#
#     def py_leeimproved(array, bounds=(0.5, 3.0), thres=5.0, looks=1.0, win=(9, 9), newsig=0.5):
#         print(bounds, thres, looks, win)
#
#         sig2 = 1.0 / looks
#         sfak = 1.0 + sig2
#         nsig2 = newsig
#         nsfak = 1.0 + nsig2
#
#         ny = array.shape[0]
#         nx = array.shape[1]
#         ym = win[0]//2
#         xm = win[1]//2
#         res = 0.0
#
#         out = np.zeros_like(array)
#         for k in range(ym, ny-ym):
#             for l in range(xm, nx-xm):
#                 m2arr = 0.0
#                 marr = 0.0
#                 n = 0
#                 for y in range(-1, 2):                          # check 3x3 neighbourhood
#                     for x in range(-1, 2):
#                         m2arr += array[k+y, l+x]**2
#                         marr += array[k+y, l+x]
#                         if array[k+y, l+x] > thres:
#                             n += 1
#                         # print(k,y, l,x, marr)
#
#                 if n >= 6:                                       # keep all point targets
#                     for y in range(-1, 2):
#                         for x in range(-1, 2):
#                             if array[k+y, l+x] > thres:
#                                 out[k+y, l+x] = array[k+y, l+x]
#
#                 if out[k, l] == 0.0:                             # no point target
#                     _m2arr = m2arr
#                     _marr = marr
#                     m2arr /= 9.0
#                     marr /= 9.0
#                     vary = (m2arr - marr**2)
#                     if vary < 1e-10: vary = 1e-10
#                     varx = ((vary - marr ** 2 * sig2) / sfak)
#                     if varx < 0: varx = 0
#                     kfac = varx / vary
#                     xtilde = (array[k, l] - marr) * kfac + marr
#                     i1 = xtilde*bounds[0]
#                     i2 = xtilde*bounds[1]
#                     m2arr = 0.0
#                     marr = 0.0
#                     n = 0
#
#                     for y in range(-ym, ym+1):
#                         for x in range(-xm, xm+1):
#                             if array[k+y, l+x]>i1 and array[k+y, l+x]<i2:
#                                 m2arr += array[k+y, l+x]**2
#                                 marr += array[k+y, l+x]
#                                 n += 1
#                     if n == 0:
#                         out[k, l] = 0.0
#                     else:
#                         m2arr /= n
#                         marr /= n
#                         vary = (m2arr - marr**2)
#                         if vary < 1e-10: vary = 1e-10
#                         varx = ((vary - marr ** 2 * nsig2) / nsfak)
#                         if varx < 0: vary = 0.0
#                         kfac = varx / vary
#                         out[k, l] = (array[k, l] - marr) * kfac + marr
#
#         return out
#
#
# except ImportError:
#     logging.info("LeeSigmaImproved cython module not found. (run build process?)")
#
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------

try:
    from .Despeckle_extensions import cy_MCB
    from .Despeckle_extensions import cy_MCBdist

    class MCB(pyrat.Worker):
        """
        Multi-Channel Beltrami Filter implementation

        :author: Joel Amao
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'Multichannel Beltrami Filter'}
        para = [
            {'var': 'orientations', 'value': 8, 'type': 'int', 'range': [2, 16],
             'text': '# of edge aligned gaussian kernels for preprocessing'},
            {'var': 'win', 'value': 7, 'type': 'int', 'range': [7, 30],
             'text': 'Window size'},
            {'var': 'looks', 'value': 4, 'type': 'float', 'range': [1, 99],
             'text': '# of looks of the input data (MLC only)'},
            {'var': 'presteps', 'value': 3, 'type': 'int', 'range': [1, 99], 'text': '# of preprocessing steps'},
            {'var': 'betastr', 'value': 2.0, 'type': 'float', 'range': [0.1, 99.0],
             'text': 'Phi (beta weight) parameter'},
            {'var': 'sigma', 'value': 1.0, 'type': 'float', 'range': [0.1, 99.0], 'text': 'Sigma parameter'},
            {'var': 'denl', 'value': 100, 'type': 'int', 'range': [1, 9999], 'text': 'Target ENL'},
            {'var': 'showlayers', 'value': False, 'type': 'bool', 'text': 'DEBUG: Shows intermediate iterations'},
            {'var': 'simdim', 'value': 100, 'type': 'int', 'range': [50, 1000], 'text': 'DEBUG: Simulated data dimensions'},
            {'var': 'dsigma', 'value': 0.0, 'type': 'float', 'range': [0.0, 99.0], 'text': 'DEBUG: Sigma change per iteration'}
        ]

        def __init__(self, *args, **kwargs):
            super(MCB, self).__init__(*args, **kwargs)
            self.name = "Multichannel Beltrami Filter"
            self.blockprocess = True
            self.blocksize = (self.win - 2) * 16
            self.blockoverlap = self.win - 2
            window_size = self.win
            neighbours_local = [1 - window_size, 1, window_size + 1, window_size, window_size - 1, - 1,
                                -window_size - 1,
                                -window_size]
            self.neighbours_local = np.asarray(neighbours_local)

        def run(self, *args, **kwargs):
            print('  Sigma:  ' + str(self.sigma) + '  Sigma increase:  ' + str(self.dsigma) + '  Phi:  ' + str(self.betastr) + '  Win:  ' + str(
                self.win) + '    # of orientations:  ' + str(self.orientations) + '   Preprocessing steps:  ' + str(
                self.presteps) + '    Looks (MLC only):  ' + str(self.looks) + '   Target ENL:   ' + str(self.denl))

            attrs = pyrat.getmeta(layer=self.layer)
            array = pyrat.getdata(layer=self.layer)
            if 'CH_pol' in attrs:  # sort a bit the channels if possible
                pol = attrs['CH_pol']
                if len(set(pol) & set(["HH", "VV", "XX"])) == 3:
                    idx = [pol.index('HH'), pol.index('VV'), pol.index('XX')]
                    array = array[idx, ...]
                    attrs['CH_pol'] = [pol[i] for i in idx]
                if len(set(pol) & set(["HH", "VV", "HV", "VH"])) == 4:
                    idx = [pol.index('HH'), pol.index('VV'), pol.index('HV'), pol.index('VH')]
                    array = array[idx, ...]
                    attrs['CH_pol'] = [pol[i] for i in idx]

            self.gsize = 2
            q = 1.5
            c = np.sqrt(self.gsize * array.shape[0] / np.pi)
            self.gdev = [c * q, c / q]

            if isinstance(array,list):
                array = np.stack(array)
            numdim = array.ndim

            if numdim > 3:
                self.preprocess = False
            else:
                self.preprocess = True

            # Simulating SAR data used to estimate the beta parameter
            sMCB, newphi = MCB.simulData(array, self.looks, self.simdim)
            if numdim == 3:  # This means the input is SLC
                sMCB = sMCB[np.newaxis, :, :, :] * sMCB[:, np.newaxis, :, :].conjugate()
            P = ProgressBar('  ' + self.name, self.denl)
            P.update(0)
            iter = 0
            enl = 1
            convergence = 200
            while enl <= self.denl and  iter < convergence:

                if enl >= numdim:
                    self.preprocess = False

                # Uses box approach for preprocessing
                self.ba = True
                print('  Current phi:' + str(self.betastr))
                if iter >= 5:# and self.preprocess == False:
                    self.sigma += self.dsigma
                sMCB, self.beta = self.simMCB(sMCB, self.sigma, self.betastr, self.win, self.preprocess, self.gdev, self.orientations,
                                              self.neighbours_local, self.ba)
                if iter != 0:
                    oldlayer = newlayer
                newlayer = self.layer_process(self.realMCB, beta=self.beta, sigma=self.sigma, betastr=self.betastr,
                                              window_size=self.win, preprocess=self.preprocess, gdev=self.gdev,
                                              orient=self.orientations, neighbours_local=self.neighbours_local, ba = self.ba)
                if self.showlayers != 1:
                    if iter != 0:
                        pyrat.delete(oldlayer, silent=True)

                pyrat.activate(newlayer, silent=True)
                enl1 = np.abs(sMCB[0, 0, ...])
                dwn = np.mean(enl1)
                up = np.sqrt(np.mean((enl1-dwn)**2))
                beta = up/dwn
                enl = 1/(beta**2)

                if enl == np.inf:
                    enl = 0
                iter += 1
                P.update(enl)
                print('  Current beta: ' + str(self.beta) + '   Current sigma:' + str(self.sigma) + '   Current ENL:  ' + str(enl))
            del P
            pyrat.activate(newlayer)
            if self.showlayers == False:
                if 'CH_pol' in attrs and (len(set(pol))<=4):
                    attrs['CH_pol'] = [p1 + p2 + '*' for p1 in attrs['CH_pol'] for p2 in attrs['CH_pol']]
                attrs['ENL'] =  enl
                pyrat.setmeta(attrs)
            else:
                pyrat.adddata(sMCB)
            pyrat.activate(newlayer, silent=True)
            return newlayer
        # Method for simulating the SAR data, the input array is used in SLC mode to calculate the oversampling weights
        # If the input is not SLC then the simulated data is multilooked using self.multilook random vectors.
        @staticmethod
        def simulData(array, looks, dim):

            numdim = array.ndim
            mu = 0
            sig = np.sqrt(.5)
            newphi = 0

            if numdim == 3:  # polarimetric vector
                array = array[..., 0:dim, 0:dim]
                maxp = dim **2
                carray = array[np.newaxis, :, :, :] * array[:, np.newaxis, :, :].conjugate()
                array =  array[np.newaxis, :, :, :]
                ldim = array.shape[0:2]
                rdim = array.shape[2:4]
                cov_array = np.rollaxis(np.rollaxis(carray, 0, start=4), 0, start=4).reshape(maxp, carray.shape[0],
                                                                                             carray.shape[1])
            else:
                array = array[..., 0:dim, 0:dim]
                ldim = array.shape[0:2]
                rdim = array.shape[2:4]
                maxp = rdim[0] * rdim[1]
                cov_array = np.rollaxis(np.rollaxis(array, 0, start=4), 0, start=4).reshape(maxp, array.shape[0],
                                                                                             array.shape[1])
            cov_mean = np.mean(cov_array, axis=0)
            cov_sqrt = sp.linalg.sqrtm(cov_mean)

            if numdim == 3:
                k1 = np.zeros((rdim[1] * rdim[0], ldim[1], ldim[0]), dtype=np.complex64)
                k2 = np.zeros((rdim[1] * rdim[0], ldim[1], ldim[0]), dtype=np.complex64)
                for i in range(maxp):
                    vectorv = np.random.normal(mu, sig, (ldim[1], 1)) + 1j * np.random.normal(mu, sig, (ldim[1], 1))
                    k1[i, :, :] = np.dot(cov_sqrt, vectorv)

                speckle_SLC = np.rollaxis(np.rollaxis(k1.reshape(rdim[0], rdim[1], ldim[0], ldim[1]), 0, start=4), 0,
                                          start=4)
                sum = np.zeros((1, 1, rdim[0], rdim[1]))
                for i in range(ldim[1]):
                    sum += np.abs(np.fft.fft2(array[0, i, ...]))
                weightSLC = sp.ndimage.filters.uniform_filter(sum, [1,1,5,5])
                OS_SLC = np.fft.ifft2(np.fft.fft2(speckle_SLC) * weightSLC)
                # Calculating the correlation coefficients for a 61x61 window
                w = 30
                ww = (2 * w + 1)
                inten = np.abs(OS_SLC[0, 0,...])
                xhat = np.mean(inten)
                m = 0
                n = 0
                e1 = np.zeros([5, 5])
                e2 = np.zeros([5, 5])
                e = np.zeros([5, 5])
                for m in range(5):
                    for n in range(5):
                        for uu in range(ww):
                            u = uu - w - 1
                            for vv in range(ww):
                                v = vv - w - 1
                                e1[m, n] += ((inten[u, v] - xhat) * (inten[u - m, v - n] - xhat))
                                e2[m, n] += (inten[u, v] - xhat) ** 2
                        e[m, n] = e1[m, n]/e2[m, n]
                print('  Correlation coefficients:')
                print(e)
                newphi = (e[0,1] + e[1,0]) / 2
                out = OS_SLC
            else:
                k1 = np.zeros((rdim[1] * rdim[0], ldim[1]), dtype=np.complex64)
                k2 = np.zeros((ldim[0], ldim[1], rdim[0], rdim[1]), dtype=np.complex64)
                res = looks - np.floor(looks)
                looks = np.int(np.floor(looks))
                for k in range(looks):
                    for i in range(maxp):
                        vectorv = np.random.normal(mu, sig, (ldim[1])) + 1j * np.random.normal(mu, sig, (ldim[1]))
                        k1[i, :] = vectorv
                    sim_SLC = np.rollaxis(k1.reshape(rdim[0], rdim[1], ldim[1]), 2, start=0)
                    k2 += sim_SLC[np.newaxis, :, :, :] * sim_SLC[:, np.newaxis, :, :].conjugate()
                k2 /= looks
                out = k2
                if res!= 0:
                    alpha = res*(0.33)/(2*looks)
                    for i in range(ldim[0]):
                        for j in range(ldim[1]):
                            mavg = out[i,j,...].reshape(rdim[0]*rdim[1])
                            mavg = np.convolve(mavg, [alpha, 1 - 2*alpha, alpha], 'same')
                            out[i, j, ...] = mavg.reshape(rdim[0],rdim[1])

            return np.squeeze(out), newphi
        @staticmethod
        def rebin(arr, *shape, **kwargs):
            """
            Imitates IDL's rebin function. Allows also phase rebining.

            :author: Andreas Reigber
            """
            phase = False  # combination of *args and fixed keywords
            if 'phase' in kwargs:
                phase = kwargs['phase']  # works only in in python 3!
            if len(shape) == 1:  # allows to pass shape as list/array or normal arguments
                if type(shape[0]) == int:
                    shape = [shape[0]]
                else:
                    shape = list(shape[0])
            oarr = arr.copy()
            oshap = oarr.shape
            for d in range(arr.ndim):
                n1 = shape[d]
                n2 = oshap[d]
                if n1 < n2:
                    s = list(oarr.shape)
                    s.insert(d + 1, n2 // n1)
                    s[d] = n1
                    if phase is True:
                        oarr = numpy.angle(numpy.exp(1j * oarr.reshape(s)).mean(d + 1))
                    else:
                        oarr = oarr.reshape(s).mean(d + 1)
                else:
                    oarr = oarr.repeat(n1 // n2, axis=d)
            return oarr
        # simMCB filters the simulated data to obtain the beta statistic used in the real MCB filter.
        @staticmethod
        def simMCB(array, sigma, betastr, window_size, preprocess, gdev, orient,  neighbours_local, ba, **kwargs):
            array = array.astype(np.complex_)
            neighM = window_size // 2 - 1
            tmpW = window_size - 2
            times = window_size // 2
            sqrt2 = np.sqrt(2)
            ldim = array.shape[0:2]
            rdim = array.shape[2:4]
            maxp = rdim[0] * rdim[1]
            cov_array = np.rollaxis(np.rollaxis(array, 0, start=4), 0, start=4).reshape(maxp, ldim[0], ldim[1])
            avg = np.zeros_like(cov_array, dtype=np.complex64)
            neighbours_global = [-rdim[1] + 1, 1, rdim[1] + 1, rdim[1]]
            all_neigh = [[0, 0]]

            #Creating the global neighbor used in region growing
            for times in range(neighM):
                for i in range(tmpW):
                    for j in range(tmpW):
                        x = window_size * (i - neighM) + (j - neighM)
                        y = rdim[1] * (i - neighM) + (j - neighM)
                        toadd = [x, y]
                        border = (times + 2)
                        if not (toadd in all_neigh) and (-border < i - neighM < border) and (
                                        -border < j - neighM < border):
                            all_neigh.append(toadd)
            all_neigh = np.asarray(all_neigh)
            distance_array, da2 = MCB.preProcessing(preprocess, 1, orient, array, cov_array, gdev, neighbours_global, ba)
            beta = np.median(da2[np.isfinite(da2)])

            # Aggressively filters the oversampled data using beta2 as threshold. The correlation coefficients are
            # still not taken into account.
            if preprocess == True:
                beta2 = np.median(distance_array[np.isfinite(distance_array)])
                distance_array[distance_array >= beta2] = np.inf

            beltramiFast = cy_MCB(cov_array, distance_array, all_neigh, neighbours_local, sigma, beta,
                                  betastr, window_size, rdim[1])
            avg = np.rollaxis(np.rollaxis(beltramiFast.reshape(rdim[0], rdim[1], ldim[0], ldim[1]), 0, start=4), 0,
                              start=4)
            return avg, beta

        # realMCB is the main part of the filter, if the input is SLC then a preprocessing step
        # is added where the data is edge filtered to better preserve statistics. This step is avoided
        # if the input is already multilooked.
        @staticmethod
        def realMCB(array, beta, sigma, betastr, window_size, preprocess, gdev, orient, neighbours_local, ba, **kwargs):
            if isinstance(array,list):
                array = np.stack(array)
            array = array.astype(np.complex_)
            numdim = array.ndim
            if numdim == 3:
                array = array[np.newaxis, :, :, :] * array[:, np.newaxis, :, :].conjugate()
            elif numdim == 2:
                array = array[np.newaxis, np.newaxis, :, :]
            array[np.isnan(array)] = 0.0

            ldim = array.shape[0:2]
            rdim = array.shape[2:4]
            maxp = rdim[0] * rdim[1]
            neighM = window_size // 2 - 1
            tmpW = window_size - 2
            times = window_size // 2

            neighbours_global = [-rdim[1] + 1, 1, rdim[1] + 1, rdim[1]]
            cov_array = np.rollaxis(np.rollaxis(array, 0, start=4), 0, start=4).reshape(maxp, ldim[0], ldim[1])
            avg = np.zeros_like(cov_array, dtype=np.complex64)
            all_neigh = [[0, 0]]
            for times in range(neighM):
                for i in range(tmpW):
                    for j in range(tmpW):
                        x = window_size * (i - neighM) + (j - neighM)
                        y = rdim[1] * (i - neighM) + (j - neighM)
                        toadd = [x, y]
                        border = (times + 2)
                        if not (toadd in all_neigh) and (-border < i - neighM < border) and (
                                        -border < j - neighM < border):
                            all_neigh.append(toadd)
            all_neigh = np.asarray(all_neigh)
            distance_array = MCB.preProcessing(preprocess, 0, orient, array, cov_array, gdev, neighbours_global, ba)
            neighbours_global2 = [-rdim[1] + 1, 1, rdim[1] + 1, rdim[1], rdim[1] - 1, -1, -rdim[1] -1, -rdim[1]]

            # Aggressively filters the oversampled data using beta2 as threshold. The correlation coefficients are
            # still not taken into account.
            if preprocess == True:
                beta2 = np.median(distance_array[np.isfinite(distance_array)])
                distance_array[distance_array >= beta2] = np.inf

            beltramiFast = cy_MCB(cov_array, distance_array, all_neigh, neighbours_local, sigma, beta,
                              betastr, window_size, rdim[1])
            avg = np.rollaxis(np.rollaxis(beltramiFast.reshape(rdim[0], rdim[1], ldim[0], ldim[1]), 0, start=4), 0,
                              start=4)
            return avg

        @staticmethod
        def matLog(matrix):
            w, v = np.linalg.eigh(matrix)
            aprime = np.dot(np.conj(v).T, np.dot(matrix, v))
            return np.dot(v, np.dot(np.diag(np.log(np.diag(aprime))), np.conj(v).T))

        @staticmethod
        def preProcessing(on, sim,  orient, array, cov_array, gdev, neighbours_global, ba):
            ldim = array.shape[0:2]
            rdim = array.shape[2:4]
            maxp = rdim[0] * rdim[1]
            orient = 2 * orient
            distance_array = np.zeros(shape=[4, maxp], dtype=np.float32) + np.inf
            distance_array2 = np.zeros(shape=[1, maxp], dtype=np.float32) + np.inf
            gdevsize = np.int(np.ceil(max(gdev)))*4 + 3
            boxapproach = ba
            #This follows the previous edge aligned kernels preprocessing approach
            if on == True and boxapproach == 0:
                gcrd = np.mgrid[:gdevsize, :gdevsize]
                gcrd = [gc - gdevsize//2 for gc in gcrd]
                karr = []
                darr = []
                ang = [(2 * np.pi / orient) * (x + 1) for x in range(orient)]
                for x in range(orient):
                    rotmat = [[np.cos(ang[x]), -np.sin(ang[x])], [np.sin(ang[x]), np.cos(ang[x])]]
                    ggcrd = [[np.dot(rotmat, [gcrd[1][x][y], gcrd[0][x][y]]) for x in range(gdevsize)] for y in range(gdevsize)]
                    ggcrd = np.rollaxis(np.asarray(ggcrd), -1)
                    ocrd = [[gc - off for gc, off in zip(ggcrd, (0, s))] for s in (-1, 0, 1)]
                    odist = [[gc / gd for gc, gd in zip(ocrd[s2], gdev)] for s2 in range(3)]
                    k = [np.exp(-0.5 * (odist[s3][0] ** 2 + odist[s3][1] ** 2)) for s3 in range(3)]
                    k = [kernel * (crd[1] >= -0.3) for kernel, crd in zip(k, ocrd)]
                    k = [kernel / kernel.sum() for kernel in k]
                    karr.append(k)
                for x in range(orient):
                    tmp = []
                    pair = (x + orient // 2) % orient
                    for y in range(3):
                        k = karr[x][y]  # y = 0 is -1 offset, 1 is central
                        if x == 2 and y == 1:
                            lastkarr = k
                        k = k[np.newaxis, np.newaxis, :]
                        numbo = 2
                        gaussie = np.copy(array)
                        data_gauss = sp.ndimage.filters.correlate(gaussie.real, weights=k) + 1j * sp.ndimage.filters.correlate(gaussie.imag, weights=k)
                        origi = np.real(gaussie[0, 0, ...])
                        if x == 2 and y == 1:
                            normal = np.real(data_gauss[0, 0, ...])
                        if x == ((2 + orient // 2) % orient) and y == numbo:
                            shifted = np.real(data_gauss[0, 0, ...])
                        dist_data = np.rollaxis(np.rollaxis(data_gauss, 0, start=4), 0, start=4).reshape(maxp,
                                                                                                         ldim[0],
                                                                                                         ldim[1])
                        tmp.append(dist_data)
                    darr.append(tmp)
                distance_par = np.zeros(shape=[maxp], dtype=np.float32)
                buff_array = np.zeros_like(cov_array)
                tempish = np.zeros_like(distance_par)
                tempish2 = np.copy(tempish)
                for k in range(maxp):
                    buff_dist = 0
                    for w in range(orient):
                        pair = (w + orient // 2) % orient
                        dist_data =  darr[w][1]
                        dist_offset = darr[pair][numbo]
                        if ldim[0] == 3:
                            distance_par = cy_MCBdist(dist_data[k, :], dist_offset[k, :])
                        else:
                            distance_par = MCB.distance_ai(dist_data[k, :], dist_offset[k, :])
                        pair_dist = distance_par
                        if pair_dist >= buff_dist:
                            buff_dist = np.copy(pair_dist)
                            tempish[k] = w
                            tempish2[k] = np.copy(pair_dist)
                            buff_array[k, :] = darr[w][1][k, :]
                if sim == 1:
                    for k in range(maxp):
                        i = 0
                        rand = np.random.randint(maxp, size=1)
                        if ldim[0] == 3:
                            distance_array2[0, k] = cy_MCBdist(buff_array[k, :], buff_array[rand[0], :])
                        else:
                            distance_array2[0, k] = MCB.distance_ai(buff_array[k, :], buff_array[rand[0], :])
                        for dx in neighbours_global:
                            nx = k + dx
                            if nx < maxp:
                                if ldim[0] == 3:
                                    distance_array[i, k] = cy_MCBdist(buff_array[k, :], buff_array[nx, :])
                                else:
                                    distance_array[i, k] = MCB.distance_ai(buff_array[k, :], buff_array[nx, :])
                                i += 1
                    # print(np.median(distance_array))


                    return np.asarray(distance_array), np.asarray(distance_array2)
                else:
                    # enl1 = np.abs(buff_array[..., 0, 0])
                    # dwn = np.mean(enl1)
                    # up = np.sqrt(np.mean((enl1 - dwn) ** 2))
                    # beta = up / dwn
                    # enl = 1 / (beta ** 2)
                    # print(enl)

                    for k in range(maxp):
                        i = 0
                        for dx in neighbours_global:
                            nx = k + dx
                            if nx < maxp:
                                if ldim[0] == 3:
                                    distance_array[i, k] = cy_MCBdist(buff_array[k, :], buff_array[nx, :])
                                else:
                                    distance_array[i, k] = MCB.distance_ai(buff_array[k, :], buff_array[nx, :])

                                i += 1
                    return np.asarray(distance_array)
            # This is the current boxcar approach
            elif on == True and boxapproach == 1:
                buff = np.copy(array)
                win = [1, 1, 3 , 3]
                out = sp.ndimage.filters.uniform_filter(buff.real, win) + 1j * sp.ndimage.filters.uniform_filter(buff.imag, win)
                buff_array = np.rollaxis(np.rollaxis(out, 0, start=4), 0, start=4).reshape(maxp,
                                                                                                         ldim[0],
                                                                                                         ldim[1])
                if sim == 1:
                    for k in range(maxp):
                        i = 0
                        rand = np.random.randint(maxp, size=1)
                        if ldim[0] == 3:
                            distance_array2[0, k] = cy_MCBdist(buff_array[k, :], buff_array[rand[0], :])
                        else:
                            distance_array2[0, k] = MCB.distance_ai(buff_array[k, :], buff_array[rand[0], :])
                        for dx in neighbours_global:
                            nx = k + dx
                            if nx < maxp:
                                if ldim[0] == 3:
                                    distance_array[i, k] = cy_MCBdist(buff_array[k, :], buff_array[nx, :])
                                else:
                                    distance_array[i, k] = MCB.distance_ai(buff_array[k, :], buff_array[nx, :])
                                i += 1
                    return np.asarray(distance_array), np.asarray(distance_array2)
                else:
                    for k in range(maxp):
                        i = 0
                        for dx in neighbours_global:
                            nx = k + dx
                            if nx < maxp:
                                if ldim[0] == 3:
                                    distance_array[i, k] = cy_MCBdist(buff_array[k, :], buff_array[nx, :])
                                else:
                                    distance_array[i, k] = MCB.distance_ai(buff_array[k, :], buff_array[nx, :])

                                i += 1
                    return np.asarray(distance_array)

            else:
                if sim == 1:
                    dist_data = cov_array
                    for k in range(maxp):
                        i = 0
                        rand = np.random.randint(maxp, size=1)
                        if ldim[0] == 3:
                            distance_array2[0, k] = cy_MCBdist(dist_data[k, :], dist_data[rand[0], :])
                        else:
                            distance_array2[0, k] = MCB.distance_ai(dist_data[k, :], dist_data[rand[0], :])
                        for dx in neighbours_global:
                            nx = k + dx
                            if nx < maxp:
                                if ldim[0] == 3:
                                    distance_array[i, k] = cy_MCBdist(dist_data[k, :], dist_data[nx, :])
                                else:
                                    distance_array[i, k] = MCB.distance_ai(dist_data[k, :], dist_data[nx, :])
                                i += 1
                    return np.asarray(distance_array), np.asarray(distance_array2)
                else:
                    dist_data = cov_array
                    for k in range(maxp):
                        i = 0
                        for dx in neighbours_global:
                            nx = k + dx
                            if nx < maxp:
                                if ldim[0] == 3:
                                    distance_array[i, k] = cy_MCBdist(dist_data[k, :], dist_data[nx, :])
                                else:
                                    distance_array[i, k] = MCB.distance_ai(dist_data[k, :], dist_data[nx, :])
                                i += 1
                    return np.asarray(distance_array)

        @staticmethod
        def distance_ai(a_mat, b_mat):
            w, v = np.linalg.eigh(a_mat) #todo: Check routines
            if np.any(w <= 0):
                return np.inf
            else:
                m1 = np.dot(v, np.dot(np.diag(np.sqrt(1. / w)), np.conj(v).T)).astype(np.complex64)
                distance = np.linalg.norm(MCB.matLog(np.dot(m1, np.dot(b_mat, m1))), 'fro')
                return distance if np.isfinite(distance) else np.inf

        @staticmethod
        def distance_le(a_mat, b_mat):
            distance = np.linalg.norm(MCB.matLog(a_mat) - MCB.matLog(b_mat), 'fro')
            return distance


    @pyrat.docstringfrom(MCB)
    def mcb(*args, **kwargs):
        return MCB(*args, **kwargs).run(*args, **kwargs)

except ImportError:
    logging.info("Multichannel Beltrami cython modules not found. (run build process?)")


try:
    from .Despeckle_extensions import cy_bilateral


    class Bilateral(pyrat.FilterWorker):
        """
        Simple bilateral speckle filter (not SAR specific). Fast implementation in Cython.

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'Bilateral (simple)'}
        para = [
            {'var': 'win', 'value': 11, 'type': 'int', 'range': [3, 999], 'text': 'Window size'},
            {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity'],
             'text': 'SAR data type'}
        ]

        def __init__(self, *args, **kwargs):
            super(Bilateral, self).__init__(*args, **kwargs)
            self.name = "BILATERAL SPECKLE FILTER"
            self.blockprocess = True
            self.blockoverlap = self.win // 2 + 1

        def filter(self, array, *args, **kwargs):
            np.seterr(divide='ignore', invalid='ignore')
            if np.iscomplexobj(array):
                array = np.abs(array)
                self.type = "amplitude"
            if self.type == "amplitude":
                array *= array
            array = cy_bilateral(array, looks=self.looks, win=self.win)
            if self.type == "amplitude":
                array = np.sqrt(array)
            array[~np.isfinite(array)] = 0.0
            np.seterr(divide='warn', invalid='warn')

            return array


    @pyrat.docstringfrom(Bilateral)
    def bilateral(*args, **kwargs):
        return Bilateral(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("Bilateral module not found. (run build process?)")


# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------


class BilateralFilter(pyrat.FilterWorker):
    """
    Bilateral speckle filter

    further information:
    D'Hondt, O., Guillaso, S., & Hellwich, O. (2013).
    Iterative bilateral filtering of polarimetric SAR data.
    Selected Topics in Applied Earth Observations and Remote Sensing, IEEE Journal of, 6(3), 1628-1639.

    Buades, A., Coll, B., & Morel, J. M. (2011).
    Non-local means denoising.
    Image Processing On Line, 2011.

    c++ impementation repository: https://github.com/odhondt/PolSAR-BLF

    :author: Toni M. del Hoyo

    # TODO:
    :param array: The image to filter (2D np.ndarray)
    :type array:
    :param win:
    :type win:
    :param looks=1.0:
    :type looks:
    :param threshold=0.5:
    :type threshold:
    :param method='original':
    :type method:
    :returns:
    """

    def __init__(self, *args, **kwargs):
        super(BilateralFilter, self).__init__(*args, **kwargs)
        self.name = "BILATERAL FILTER"
        self.blockprocess = True
        self.blocksize = 30
        if 'gammaS' not in self.__dict__:
            self.gammaS = 2.8
        if 'gammaR' not in self.__dict__:
            self.gammaR = 1.33
        if 'nit' not in self.__dict__:
            self.nit = 1
        self.H = np.int(np.ceil(np.sqrt(3) * self.gammaS))
        self.blockoverlap = self.H

    def filter(self, array, *args, **kwargs):
        n, n, nry, nrx = array.shape
        out = np.zeros(array.shape, dtype=np.complex64)
        gr2 = self.gammaR * self.gammaR
        gs2 = self.gammaS * self.gammaS

        # GENERATE GAUSSIAN WINDOW
        win = np.zeros((2 * self.H + 1, 2 * self.H + 1), dtype=float)
        for x in range(0, 2 * self.H + 1):
            for y in range(0, 2 * self.H + 1):
                t = y - self.H
                s = x - self.H
                win[y, x] = np.exp(-(s ** 2 + t ** 2) / gs2)

        for it in range(self.nit):
            # FOR EVERY PIXEL IN THE IMAGE
            for xpos in range(nrx):
                for ypos in range(nry):

                    #  FOR EVERY PIXEL IN THE WINDOW (EXCEPT CENTRAL)
                    tmin = np.maximum(ypos - self.H, 0)
                    smin = np.maximum(xpos - self.H, 0)
                    tmax = np.minimum(ypos + self.H, nry)
                    smax = np.minimum(xpos + self.H, nrx)

                    mat1 = np.linalg.inv(sp.linalg.sqrtm(array[..., ypos, xpos]))
                    wr_max = 0
                    wsum = 0
                    filt = np.zeros((n, n), dtype=np.complex64)

                    for s in range(smin, smax):
                        for t in range(tmin, tmax):
                            wr = 0
                            # COMPUTE DISTANCE D BW. CENTRAL PIXEL AND OTHERS
                            d = np.linalg.norm(np.dot(np.dot(mat1, array[..., t, s]), mat1), ord='fro')
                            # d = np.linalg.norm(array[..., t, s] - array[..., ypos, xpos], ord='fro')

                            # TRANSLATE D INTO RADIOMETRIC WEIGHT WR
                            if not np.isnan(d):
                                wr = np.exp(-d / gr2)

                            # UPDATE MAXIMUM WEIGHT WRMAX
                            if wr_max < wr < 1:
                                wr_max = wr

                            # GENERATE SPATIAL DISTANCE WEIGHT WD
                            wd = win[t - ypos + self.H, s - xpos + self.H]

                            # GENERATE GLOBAL WEIGHT W = WR * WD
                            w = wr * wd

                            # FILTER PIXEL
                            filt += w * array[..., t, s]

                            # UPDATE SUM OF WEIGHTS WSUM
                            wsum += w

                    # FILTER CENTRAL PIXEL BY WRMAX
                    filt += wr_max * array[..., ypos, xpos]

                    # UPDATE WSUM
                    wsum += wr_max

                    # DIVIDE FILTERED PIXEL BY WSUM
                    if wsum > 1e-10:
                        filt /= wsum
                    else:
                        filt = array[..., ypos, xpos]

                    out[..., ypos, xpos] = filt

            # UPDATE INPUT FOR NEXT ITERATION
            array = out
        return out


@pyrat.docstringfrom(BilateralFilter)
def bilateralfilter(*args, **kwargs):
    return BilateralFilter(*args, **kwargs).run(**kwargs)


# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------

try:
    from pyrat.tools import ProgressBar
    from .Despeckle_extensions import cy_srad


    class SRAD(pyrat.Worker):
        """
        Anisotropic diffusion speckle filter. Fast implementation in Cython.
        Y. You and S.T. Acton: Speckle Reducing Anisotropic diffusion
        IEEE Transactions on Image Processing, Vol. 11, No.11, pp. 1260-1269-213, 2002']

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'Anisotropic diffusion'}
        para = [
            {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'step', 'value': 0.15, 'type': 'float', 'range': [0.0, 1.0], 'text': 'step size'},
            {'var': 'iter', 'value': 100, 'type': 'int', 'range': [2, 1000], 'text': '# of iterations'},
            {'var': 'scale', 'value': 1.0, 'type': 'float', 'range': [1.0, 10.0], 'text': 'scale'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity'],
             'text': 'SAR data type'}]

        def __init__(self, *args, **kwargs):
            super(SRAD, self).__init__(*args, **kwargs)
            self.name = "SRAD SPECKLE FILTER"
            self.blockprocess = True
            self.blockoverlap = 1

        def run(self, *args, **kwargs):
            P = ProgressBar('  ' + self.name, self.iter)
            P.update(0)
            for k in range(self.iter):
                if k != 0:
                    oldlayer = newlayer
                newlayer = self.layer_process(self.srad, looks=self.looks, step=self.step, iter=k,
                                              scale=self.scale, type=self.type)
                if k != 0:
                    pyrat.delete(oldlayer, silent=True)
                pyrat.activate(newlayer, silent=True)
                P.update(k + 1)
            del P
            pyrat.activate(newlayer)
            return newlayer

        @staticmethod
        def srad(array, looks=1, step=0.05, iter=0, scale=1.0, type='amplitude', **kwargs):
            if array.ndim == 3:  # polarimetric vector
                if np.iscomplexobj(array) or type == "amplitude":
                    array = np.abs(array) ** 2
                    type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                type = "amplitude"
            elif array.ndim == 4:  # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                type = "intensity"
            else:  # single channel data
                if np.iscomplexobj(array) or type == "amplitude":
                    array = np.abs(array) ** 2
                    type = "amplitude"
                else:
                    array = np.abs(array)
                span = array.copy()
                array = array[np.newaxis, np.newaxis, ...]
            array = cy_srad(span.clip(1e-10), array.clip(1e-10), looks=looks, step=step, iter=iter, scale=scale)
            array[~np.isfinite(array)] = 0.0
            if type == "amplitude":
                array[array < 0] = 0.0
                array = np.sqrt(array)
            return np.squeeze(array)


    @pyrat.docstringfrom(SRAD)
    def srad(*args, **kwargs):
        return SRAD(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("SRAD cython module not found. (run build process?)")

# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------
# try:
#     from pyrat.tools import ProgressBar
#     from .Despeckle_extensions import cy_srad_old
#
#     class SRAD_OLD(pyrat.Worker):
#         """
#         Anisotropic diffusion speckle filter...
#         Y. You and S.T. Acton: Speckle Reducing Anisotropic diffusion
#         IEEE Transactions on Image Processing, Vol. 11, No.11, pp. 1260-1269-213, 2002']
#
#         :author: Andreas Reigber
#         """
#         gui = {'menu': 'SAR|Speckle filter', 'entry': 'Diffusion'}
#         para = [
#             {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
#             {'var': 'step', 'value': 0.15, 'type': 'float', 'range': [0.0, 1.0], 'text': 'step size'},
#             {'var': 'iter', 'value': 100, 'type': 'int', 'range': [2, 1000], 'text': '# of iterations'},
#             {'var': 'scale', 'value': 1.0, 'type': 'float', 'range': [1.0, 10.0], 'text': 'scale'}            ]
#
#         def __init__(self, *args, **kwargs):
#             super(SRAD_OLD, self).__init__(*args, **kwargs)
#             self.name = "SRAD SPECKLE FILTER"
#             self.blockprocess = True
#             self.blockoverlap = 1
#             # self.nthreads = 1
#
#         def run(self, *args, **kwargs):
#             P = ProgressBar('  ' + self.name, self.iter)
#             P.update(0)
#
#             look = [self.looks]
#
#             for k in range(self.iter):
#                 if k != 0:
#                     oldlayer = newlayer
#                 newlayer = self.layer_process(self.srad, looks=self.looks, step=self.step, iter=k, scale=self.scale)
#                 if k != 0:
#                     pyrat.delete(oldlayer, silent=True)
#                 pyrat.activate(newlayer, silent=True)
#
#                 # data = pyrat.getdata()[50:150, 550:650]
#                 # data = pyrat.getdata()[25:75, 275:325]
#                 # look.append(np.mean(data)**2 / np.var(data))
#                 # print(look[-1])
#                 # self.looks = look[-1]
#
#                 P.update(k+1)
#             del P
#
#             # q0 = [np.sqrt(look[0]/l) for l in look]
#             # q1 = np.exp(-self.step*np.arange(self.iter)/6.0)
#             # # q0the = 1.0/np.sqrt(self.looks)*np.exp(-self.step*np.arange(self.iter)/4.0)
#             #
#             # # plot(q0)
#             # # plot(q1)
#             # stop()
#             pyrat.activate(newlayer)
#             return newlayer
#
#         @staticmethod
#         def srad(array, looks=1, step=0.05, iter=0, scale=1.0, **kwargs):
#             out = cy_srad_old(array.clip(1e-10), looks=looks, step=step, iter=iter, scale=scale)
#             return out
#
#     def srad_old(*args, **kwargs):
#         return SRAD_OLD(*args, **kwargs).run(**kwargs)
#
# except ImportError:
#     logging.info("SRAD cython module not found. (run build process?)")
#
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------

try:
    from .Despeckle_extensions import cy_emdes as cy_emdes


    class EMDES(pyrat.FilterWorker):
        """
        Test filter

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'EMDES'}
        para = [
            {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
             'subtext': ['range', 'azimuth']},
            {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity'],
             'text': 'SAR data type'}
        ]

        def __init__(self, *args, **kwargs):
            super(EMDES, self).__init__(*args, **kwargs)
            self.name = "EMDES SPECKLE FILTER"
            self.blockprocess = True
            self.blockoverlap = self.win[0] // 2 + 1

        def filter(self, array, *args, **kwargs):
            if array.ndim == 3:  # polarimetric vector
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array) ** 2
                    self.type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                self.type = "amplitude"
            elif array.ndim == 4:  # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                self.type = "intensity"
            else:  # single channel data
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array) ** 2
                    self.type = "amplitude"
                else:
                    array = np.abs(array)
                span = array.copy()
                array = array[np.newaxis, np.newaxis, ...]
            array = cy_emdes(span, array, looks=self.looks, win=self.win)
            array[~np.isfinite(array)] = 0.0
            if self.type == "amplitude":
                array[array < 0] = 0.0
                array = np.sqrt(array)
            return np.squeeze(array)


    @pyrat.docstringfrom(EMDES)
    def emdes(*args, **kwargs):
        return EMDES(*args, **kwargs).run(**kwargs)

except ImportError:
    logging.info("EMDES cython module not found. (run build process?)")

try:
    from .Despeckle_extensions import cy_idanq as cy_idanq
    from scipy.ndimage.filters import median_filter


    class IDANQ(pyrat.FilterWorker):
        """
        IDAN Speckle Filter (Intensity Driven Adaptive Neighbourhood).

        This is a "quick'n'dirty" implementation, using only the span
        of multidimensional data for region growing. Fast implementation in Cython.

        for details see: G. Vasile et al.: "Intensity-Driven Adaptive-Neighborhood
        Technique for Polarimetric and Interferometric SAR Parameters Estimation",
        IEEE Trans. Geosc. Remote Sensing, 44(6), 2006

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'IDAN quick'}
        para = [
            {'var': 'size', 'value': 50, 'type': 'int', 'range': [3, 999], 'text': 'Neighbourhood size'},
            {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity'],
             'text': 'SAR data type'},
            {'var': 'llmmse', 'value': True, 'type': 'bool', 'text': 'LLMMSE filtering'}
        ]

        def __init__(self, *args, **kwargs):
            super(IDANQ, self).__init__(*args, **kwargs)
            self.name = "IDAN-Q SPECKLE FILTER"
            self.blockprocess = True
            # self.nthreads = 1
            self.blockoverlap = self.size + 1
            self.blocksize = self.blockoverlap * 4

        def filter(self, array, *args, **kwargs):
            if array.ndim == 3:  # polarimetric vector
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array) ** 2
                    self.type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                self.type = "amplitude"
            elif array.ndim == 4:  # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                self.type = "intensity"
            else:  # single channel data
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array) ** 2
                    self.type = "amplitude"
                else:
                    array = np.abs(array)
                span = array.copy()
                array = array[np.newaxis, np.newaxis, ...]
            array = cy_idanq(span, array, looks=self.looks, nmax=self.size, llmmse=self.llmmse)
            array[~np.isfinite(array)] = 0.0
            if self.type == "amplitude":
                array[array < 0] = 0.0
                array = np.sqrt(array)
            return np.squeeze(array)


    @pyrat.docstringfrom(IDANQ)
    def idanq(*args, **kwargs):
        return IDANQ(*args, **kwargs).run(*args, **kwargs)

except ImportError:
    logging.info("IDANQ cython module not found. (run build process?)")

