import pyrat
import scipy as sp
from scipy import optimize as opt
import numpy as np
import logging
from scipy.ndimage import filters


class Boxcar(pyrat.FilterWorker):
    """
    Boxcar / Moving average (speckle) filter.

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


class Gauss(pyrat.FilterWorker):
    """
    Gaussian (speckle) filter...

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


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

try:
    from .Despeckle_extensions import cy_leesigma

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
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity'], 'text': 'SAR data type'}
            ]

        def __init__(self, *args, **kwargs):
            super(LeeSigma, self).__init__(*args, **kwargs)
            self.name = "SIGMA LEE SPECKLE FILTER"
            self.blockprocess = True
            self.blockoverlap = self.win[0] // 2 + 1

        def filter(self, array, *args, **kwargs):
            if array.ndim == 3:                                                # polarimetric vector
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array)**2
                    self.type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                self.type = "amplitude"
            elif array.ndim == 4:                                              # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                self.type = "intensity"
            else:                                                              # single channel data
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array)**2
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
        Lee's improved sigma speckle filter...

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
            bounds = opt.fmin(self.optf, [0.5, 2.0], args=(self.looks, self.sigma), disp=False)    # calc sigma bounds
            newsig = self.newsig(bounds[0], bounds[1], sigrng=self.sigma, looks=self.looks)        # calc new stddev
            P.update(0)
            perc = 100.0-self.perc*100.0                                                       # point target theshold
            pthreshold = np.mean(self.layer_accumulate(self.estimate_percentile, type=self.type, perc=perc))
            P.update(2)

            layer = self.layer_process(self.leeimproved, bounds=bounds, newsig=newsig,
                                       thres=pthreshold, looks=self.looks, win=self.win, type=self.type)
            P.update(10)
            del P
            pyrat.activate(layer)
            return layer

        @staticmethod
        def leeimproved(array, bounds=(1, 2), thres=10.0, looks=1.0, win=(7, 7), newsig=0.5, type='amplitude', **kwargs):

            if array.ndim == 3:                                                # polarimetric vector
                if np.iscomplexobj(array) or type == "amplitude":
                    array = np.abs(array)**2
                    type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                type = "amplitude"
            elif array.ndim == 4:                                              # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                type = "intensity"
            else:                                                              # single channel data
                if np.iscomplexobj(array) or type == "amplitude":
                    array = np.abs(array)**2
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
            if array.ndim == 3:                                                # polarimetric vector
                if np.iscomplexobj(array) or type == "amplitude":
                    span = np.sum(np.abs(array)**2, axis=0)
                else:
                    span = np.sum(np.abs(array), axis=0)
            elif array.ndim == 4:                                              # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
            else:                                                              # single channel data
                if np.iscomplexobj(array) or type == "amplitude":
                    span = np.abs(array)**2
                else:
                    span = np.abs(array)
            return np.percentile(span, perc)

        @staticmethod
        def specklepdf(i, looks=1.0):
            if i < 0.0:
                return 0.0
            else:
                return ((looks**looks) * (i**(looks-1.0))) / sp.special.gamma(looks) * np.exp(-looks*i)

        @staticmethod
        def meanpdf(i, looks=1.0):
            if i < 0.0:
                return 0.0
            else:
                return ((looks**looks) * (i**(looks-1.0))) / sp.special.gamma(looks) * np.exp(-looks*i) * i

        def sigpdf(self, i, looks=1.0):
            return (i - 1.0)**2 * self.specklepdf(i, looks=1.0)

        def newsig(self, i1, i2, sigrng=0.9, looks=1.0):
            return 1 / sigrng * sp.integrate.quad(self.sigpdf, i1, i2, args=(looks, ))[0]

        def sigmarange(self, i1, i2, looks=1.0):
            return np.clip(sp.integrate.quad(self.specklepdf, i1, i2, args=(looks, ))[0], 1e-10, 1.0)

        def intmean(self, i1, i2, looks=1.0):
            return 1.0 / self.sigmarange(i1, i2, looks) * sp.integrate.quad(self.meanpdf, i1, i2, args=(looks, ))[0]

        def optf(self, i, looks, sigr):
            return (self.sigmarange(i[0], i[1], looks) - sigr)**2 + (self.intmean(i[0], i[1], looks) - 1.0)**2


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
    from .Despeckle_extensions import cy_bilateral

    class Bilateral(pyrat.FilterWorker):
        """
        Simple bilateral speckle filter

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'Bilateral (simple)'}
        para = [
            {'var': 'win', 'value': 11, 'type': 'int', 'range': [3, 999], 'text': 'Window size'},
            {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude','intensity'], 'text': 'SAR data type'}
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
        Anisotropic diffusion speckle filter...
        Y. You and S.T. Acton: Speckle Reducing Anisotropic diffusion
        IEEE Transactions on Image Processing, Vol. 11, No.11, pp. 1260-1269-213, 2002']

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'SRAD'}
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
            self.nthreads = 1

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
                P.update(k+1)
            del P
            pyrat.activate(newlayer)
            return newlayer

        @staticmethod
        def srad(array, looks=1, step=0.05, iter=0, scale=1.0, type='amplitude', **kwargs):
            if array.ndim == 3:                                                # polarimetric vector
                if np.iscomplexobj(array) or type == "amplitude":
                    array = np.abs(array)**2
                    type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                type = "amplitude"
            elif array.ndim == 4:                                              # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                type = "intensity"
            else:                                                              # single channel data
                if np.iscomplexobj(array) or type == "amplitude":
                    array = np.abs(array)**2
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
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

try:
    from .Despeckle_extensions import cy_emdes as cy_emdes
    class EMDES(pyrat.FilterWorker):
        """
        Test filter

        :author: Andreas Reigber
        """
        gui = {'menu': 'SAR|Speckle filter', 'entry': 'EMDES'}
        para = [
            {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size', 'subtext': ['range', 'azimuth']},
            {'var': 'looks', 'value': 2.0, 'type': 'float', 'range': [1.0, 99.0], 'text': '# of looks'},
            {'var': 'type', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity'], 'text': 'SAR data type'}
            ]

        def __init__(self, *args, **kwargs):
            super(EMDES, self).__init__(*args, **kwargs)
            self.name = "EMDES SPECKLE FILTER"
            self.blockprocess = True
            self.blockoverlap = self.win[0] // 2 + 1

        def filter(self, array, *args, **kwargs):
            if array.ndim == 3:                                                # polarimetric vector
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array)**2
                    self.type = "amplitude"
                else:
                    array = np.abs(array)
                span = np.sum(array, axis=0)
                array = array[np.newaxis, ...]
                self.type = "amplitude"
            elif array.ndim == 4:                                              # covariance data
                span = np.abs(np.trace(array, axis1=0, axis2=1))
                self.type = "intensity"
            else:                                                              # single channel data
                if np.iscomplexobj(array) or self.type == "amplitude":
                    array = np.abs(array)**2
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

