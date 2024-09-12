import matplotlib.pyplot as plt

import numpy as np
from scipy.ndimage import filters
from functools import reduce

import pyrat
from pyrat.lib.ste import Blocxy
from pyrat.tools import linspace_nd
import scipy.signal as ss

def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

class Boxcar(pyrat.FilterWorker):
    """
    Boxcar / Moving average phase noise filter.

    :author: Andreas Reigber
    """

    gui = {'menu': 'InSAR|Phase noise filter', 'entry': 'Boxcar'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
         'subtext': ['range', 'azimuth']},
    ]

    def __init__(self, *args, **kwargs):
        super(Boxcar, self).__init__(*args, **kwargs)
        self.name = "BOXCAR FILTER"
        self.blockprocess = True
        self.blockoverlap = self.win[0] // 2 + 1
        self.scaling_hint = 'phase'

    def filter(self, array, *args, **kwargs):
        win = self.win
        if array.ndim == 3:
            win = [1] + self.win
        if array.ndim == 4:
            win = [1, 1] + self.win
        array[np.isnan(array)] = 0.0
        if np.iscomplexobj(array):
            return filters.uniform_filter(array.real, win) + 1j * filters.uniform_filter(array.imag, win)
        else:
            tmp = np.exp(1j * array)
            tmp = filters.uniform_filter(tmp.real, win) + 1j * filters.uniform_filter(tmp.imag, win)
            return np.angle(tmp)


@pyrat.docstringfrom(Boxcar)
def boxcar(*args, **kwargs):
    return Boxcar(*args, **kwargs).run(*args, **kwargs)


class Goldstein(pyrat.FilterWorker):
    """
    Goldstein phase noise filter

    further information:
    R.M. Goldstein and C.L. Werner: Radar interferogarm filtering
    for geophysical applications, Geophys. Res. Lett., Vol. 25,
    No. 21, pp. 4035-4038, 1998

    :author: Andreas Reigber
    """
    gui = {'menu': 'InSAR|Phase noise filter', 'entry': 'Goldstein'}
    para = [
        {'var': 'exp', 'value': 0.5, 'type': 'float', 'text': 'Spectral exponent'},
        {'var': 'bs', 'value': 32, 'type': 'int', 'text': 'Window size'},
        {'var': 'smo', 'value': False, 'type': 'bool', 'text': 'Perform a prior gaussian smoothing'}
    ]

    def __init__(self, *args, **kwargs):
        super(Goldstein, self).__init__(*args, **kwargs)
        self.name = "GOLDSTEIN FILTER"
        self.blockprocess = True
        self.blockoverlap = self.bs
        self.scaling_hint = 'phase'

    def filter(self, array, *args, **kwargs):
        if self.win is None:
            self.win = [3, 1]

        array = np.exp(1j * array)
        bp = Blocxy(array, (self.bs, self.bs), (self.bs // 2, self.bs // 2))
        for block in bp.getiterblocks():
            block = np.fft.fft2(block)
            if not self.smo:
                block *= abs(block) ** self.exp
            else:
                block *= abs(self.smooth(block, win=self.win)) ** self.exp
            block = np.fft.ifft2(block)
            bp.setiterblocks(block)
        output = bp.getresult()
        output = np.angle(output)
        return output

    @staticmethod
    def smooth(array, win):
        return filters.uniform_filter(array.real, win) + 1j * filters.uniform_filter(array.imag, win)

@pyrat.docstringfrom(Goldstein)
def goldstein(*args, **kwargs):
    return Goldstein(*args, **kwargs).run(*args, **kwargs)

class Goldstein_iter(pyrat.FilterWorker):
    """
    An iterative Goldstein SAR interferogram filter

    further information:
    Zhao, Chaoying & Zhang, Qin & Ding, Xiaoli & Zhang, Jing. (2012).
    An iterative Goldstein SAR interferogram filter. International Journal
    of Remote Sensing. 33. 3443-3455. 10.1080/01431161.2010.532171.
    "https://www.researchgate.net/publication/233269737_An_iterative_Goldstein_SAR_interferogram_filter"

    :arg array: The image to filter (numpy.ndarray)
    :param exp: forces a filter strength. The given value will force a filtering strength. Defaults to None.
    :param bsx: window size in first dimension. Defaults to 32.
    :param bsy: window size in second dimension. Defaults to 32.
    :param prc: patch overlap (percentual). Defaults to 10%.
    :param tol: filtering tolerance for stopping the filter. Defaults to 0.1.
    :param smo: performs a smoothing on the spectra before filtering. Defaults to off.
    :param win: smoothing window size. Defaults to [3,1].
    :param force_coh: forces a (1-coh) value instead of the exponential value. Off by default.
    :param force_coh2: forces a (1-coh**2) value instead of the exponential value. Off by default.
    :param no_iter: do not iterate. Off by default.

    :author: Joel Amao, Andreas Reigber, Muriel Pinheiro
    """
    gui = {'menu': 'InSAR|Phase noise filter', 'entry': 'Goldstein'}
    para = [
        {'var': 'exp', 'value': None, 'type': 'float', 'text': 'Spectral exponent'},
        {'var': 'bsx', 'value': 32, 'type': 'int', 'text': 'Window size (x dim)'},
        {'var': 'bsy', 'value': 32, 'type': 'int', 'text': 'Window size (y dim)'},
        {'var': 'prc', 'value': 0.5, 'type': 'float', 'text': 'Patch overlap (%)'},
        {'var': 'tol', 'value': 0.1, 'type': 'float', 'text': 'Filtering tolerance'},
        {'var': 'smo', 'value': None, 'type': 'bool', 'text': 'Perform a prior gaussian smoothing'},
        {'var': 'force_coh', 'value': None, 'type': 'bool', 'text': 'forces a (1-coh) value instead of the exponential '
                                                                    'value. Off by default'},
        {'var': 'force_coh2', 'value': None, 'type': 'bool', 'text': 'forces a (1-coh**2) value instead of the '
                                                                     'exponential value. Off by default'},
        {'var': 'no_iter', 'value': None, 'type': 'bool', 'text': 'Turns off iterations. Off by default'}
    ]


    def __init__(self, *args, **kwargs):
        super(Goldstein_iter, self).__init__(*args, **kwargs)
        self.name = "ITERATIVE GOLDSTEIN FILTER"
        self.blockprocess = True
        self.blockoverlap = self.bsx
        self.scaling_hint = 'phase'
        self.win = None

    def filter(self, array, *args, **kwargs):

        if self.win is None:
            self.win = [3, 1]

        array = np.exp(1j * array)
        bp = Blocxy(array, (self.bsx, self.bsy), (self.bsx // int((1 // self.prc)), self.bsy // int((1 // self.prc))))

        for block in bp.getiterblocks():
            pc = (np.abs(np.sum(block))) / (np.sum(np.abs(block)))
            alpha = 1 - pc
            pc_old = pc - 2 * self.tol

            if self.exp is not None:
                alpha = self.exp
            elif self.force_coh is not None:
                coh = np.abs(self.smooth(np.exp(-1j * block), [3, 3]))
                coh = np.clip(np.nan_to_num(coh), 0.0, 1.0)  # get rid of numerical inaccuracies!
                alpha = 1 - np.mean(coh)
            elif self.force_coh2 is not None:
                coh = np.abs(self.smooth(np.exp(-1j * block), [3, 3]))
                coh = np.clip(np.nan_to_num(coh), 0.0, 1.0)  # get rid of numerical inaccuracies!
                alpha = 1 - np.mean(coh)**2

            if self.exp is not None or self.force_coh is not None or self.force_coh2 is not None:
                self.no_iter = True

            while pc - pc_old >= self.tol:
                block = np.fft.fft2(block)
                if self.smo:
                    block *= abs(self.smooth(block, win=self.win)) ** alpha
                else:
                    block *= abs(block) ** alpha
                block = np.fft.ifft2(block)
                pc_old = pc
                pc = (np.abs(np.sum(block))) / (np.sum(np.abs(block)))
                if self.no_iter is None:
                    if pc >= pc_old:
                        bp.setiterblocks(block)
                        alpha = 1 - pc
                else:
                    bp.setiterblocks(block)
                    pc = pc_old

        output = bp.getresult()
        output = np.angle(output)
        return output

    @staticmethod
    def smooth(array, win):
        return filters.uniform_filter(array.real, win) + 1j * filters.uniform_filter(array.imag, win)

@pyrat.docstringfrom(Goldstein_iter)
def goldstein_iter(*args, **kwargs):
    return Goldstein_iter(*args, **kwargs).run(*args, **kwargs)

class Goldstein_lfe(pyrat.FilterWorker):
    """
    Improved Goldstein Interferogram Filter Based on Local Fringe Frequency Estimation

    further information:
    Qingqing Feng, Huaping Xu, Zhefeng Wu, Yanan You, Wei Liuand Shiqi Ge (2016)
    Improved Goldstein Interferogram Filter Based on Local Fringe Frequency Estimation
    Sensors 2016, 16(11), 1976; https://doi.org/10.3390/s16111976
    "https://www.mdpi.com/1424-8220/16/11/1976"

    :author: Joel Amao
    """
    gui = {'menu': 'InSAR|Phase noise filter', 'entry': 'Goldstein LFE'}
    para = [
        {'var': 'exp', 'value': None, 'type': 'float', 'text': 'Spectral exponent'},
        {'var': 'bs', 'value': 32, 'type': 'int', 'text': 'Window size'},
        {'var': 'prc', 'value': 0.5, 'type': 'float', 'text': 'Patch overlap (%)'},
        {'var': 'force_coh', 'value': True, 'type': 'bool', 'text': 'forces a (1-coh) value instead of the exponential '
                                                                    'value. On by default'},
        {'var': 'force_coh2', 'value': False, 'type': 'bool', 'text': 'forces a (1-coh**2) value instead of the '
                                                                     'exponential value. Off by default'},
        {'var': 'win', 'value': [3, 3], 'type': 'int', 'range': [3, 999], 'text': 'Window size'}
    ]


    def __init__(self, *args, **kwargs):
        super(Goldstein_lfe, self).__init__(*args, **kwargs)
        self.name = "LOCAL FRINGE GOLDSTEIN FILTER "
        self.blockprocess = True
        self.blockoverlap = self.bs
        self.scaling_hint = 'phase'
        self.bsx = self.bsy = self.bs
        self.text = 'None'
        self.debug = False
    def filter(self, array, *args, **kwargs):

        def get_res(block, pha):
            # This function takes the smoothed spectra and an interferogram to estimate  local slopes
            # Returns row and col freqs. and max autocorr value

            corr = np.abs(block)
            size = corr.shape[0]
            amax = np.unravel_index(corr.argmax(), corr.shape)
            row = np.fft.fftfreq(size)[amax[0]]
            col = np.fft.fftfreq(size)[amax[1]]
            return row, col, np.max(corr)

        if self.win is None:
            self.win = [3, 3]

        if self.force_coh2:
            self.force_coh = False

        iscomplex = True
        if not np.iscomplexobj(array):
            # Getting interferogram from phase array
            iscomplex = False
            array = np.exp(1j * array)

        bp = Blocxy(array, (self.bsx, self.bsy), (self.bsx // int((1 // self.prc)), self.bsy // int((1 // self.prc))))

        for block in bp.getiterblocks():
            # Storing original phase interferogram
            org_pha = block.copy()
            block = np.exp(1j * np.angle(block))

            # Interferogram spectrum and smoothing
            block = np.fft.fft2(block)
            block = self.smooth(block, win=self.win)

            # Gets residual interferogram
            freqx, freqy, peak_val = get_res(block, org_pha)

            # Removes estimated fringes from original interferogram
            x = np.linspace(0, self.bs - 1, self.bs) * freqx
            y = np.linspace(0, self.bs - 1, self.bs) * freqy
            xv, yv = np.meshgrid(y, x)
            res_pha = org_pha * np.exp(-1j * 2 * np.pi * (xv + yv))
            # res_pha = [[org_pha[i, j] * np.exp(-1j * 2 * np.pi * (freqx * i + freqy * j)) for j in range(self.bs)] for i in
            #        range(self.bs)]
            # res_pha = np.asarray(res_pha)

            res_spc = np.fft.fft2(res_pha)
            coh = peak_val/(self.bs**2)

            freqx_res, freqy_res, peak_val = get_res(res_spc, res_pha)

            if self.exp is not None:
                alpha = self.exp
            elif self.force_coh and self.exp is None:
                alpha = 1 - coh
            elif self.force_coh2 and self.exp is None:
                alpha = 1 - coh**2

            alpha += np.sqrt(np.abs(freqx_res)**2 + np.abs(freqy_res)**2)

            # Filtering residual spectrum
            res_spc *= np.abs(res_spc) ** alpha
            fil_pha = np.fft.ifft2(res_spc)

            # Adding back fringes to the complex interferogram
            fil_pha = [[fil_pha[i, j] * np.exp(1j * 2 * np.pi * (freqx * i + freqy * j)) for j in range(self.bs)]
                        for i in range(self.bs)]
            fil_pha = np.asarray(fil_pha)
            bp.setiterblocks(fil_pha)
        output = bp.getresult()
        if not iscomplex:
            output = np.angle(output)

        return output

    @staticmethod
    def smooth(array, win):
        return filters.uniform_filter(array.real, win) + 1j * filters.uniform_filter(array.imag, win)

@pyrat.docstringfrom(Goldstein_lfe)
def goldstein_lfe(*args, **kwargs):
    return Goldstein_lfe(*args, **kwargs).run(*args, **kwargs)

class LFFE_boxcar(pyrat.FilterWorker):
    """
    (Local Fringe Frequency Estimator) Boxcar filter with local fringe removal

    :author: Joel Amao
    """
    gui = {'menu': 'InSAR|Phase noise filter', 'entry': 'Local Fringe Boxcar'}
    para = [
        {'var': 'bs', 'value': 32, 'type': 'int', 'text': 'Window size'},
        {'var': 'prc', 'value': 0.5, 'type': 'float', 'text': 'Patch overlap (%)'},
        {'var': 'win', 'value': [3, 3], 'type': 'int', 'range': [3, 999], 'text': 'Window size'},
        {'var': 'awin', 'value': False, 'type': 'bool', 'text': 'Uses adaptive windows based on coherence'},
        {'var': 'debug', 'value': False, 'type': 'bool', 'text': 'for debugging'},

    ]

    def __init__(self, *args, **kwargs):
        super(LFFE_boxcar, self).__init__(*args, **kwargs)
        self.name = "LFFE BOXCAR PHASE FILTER "
        self.blockprocess = True
        self.blockoverlap = self.bs
        self.scaling_hint = 'phase'
        self.bsx = self.bsy = self.bs
        self.text = 'None'

    def filter(self, array, *args, **kwargs):

        def get_res(block, pha, bs):
            # Gets residual (flat) phase and autocorrelation peak
            corr = np.abs(block)
            size = corr.shape[0]
            amax = np.unravel_index(corr.argmax(), corr.shape)
            row = np.fft.fftfreq(size)[amax[0]]
            col = np.fft.fftfreq(size)[amax[1]]
            res = np.zeros_like(pha)

            # b defines the 3-dB value
            coh = np.max(corr)/(bs**2)
            if coh < 0.25:
                # print("Coherence in patch " + str(coh))
                row = 0
                col = 0

            res = [[pha[i,j] * np.exp(-1j * 2 * np.pi * (row*i + col*j)) for j in range(size)] for i in range(size)]
            res = np.asarray(res)
            return res, row, col, np.max(corr)

        if self.win is None:
            self.win = [3, 3]

        iscomplex = True
        if not np.iscomplexobj(array):
            # Getting interferogram from phase array
            iscomplex = False
            array = np.exp(1j * array)

        bp = Blocxy(array, (self.bsx, self.bsy), (self.bsx // int((1 // self.prc)), self.bsy // int((1 // self.prc))))

        for block in bp.getiterblocks():
            org_pha = block.copy()
            # org_pha = np.exp(1j * np.angle(block))

            # Interferogram spectrum
            block = np.exp(1j * np.angle(block))
            block = np.fft.fft2(block)

            # Gets residual interferogram
            res_pha, freqx, freqy, peak_val = get_res(block, org_pha, self.bs)

            if self.awin:
                coh = peak_val / (self.bs ** 2)
                win = [min(int(1/coh), int(self.win[0]*2)), min(int(1/coh), int(self.win[1]*2))]
            else:
                win = self.win

            # Phase smoothing
            fil_pha = self.smooth(res_pha, win=win)

            # Adds fringes back
            end_pha = [[fil_pha[i, j] * np.exp(1j * 2 * np.pi * (freqx * i + freqy * j)) for j in range(self.bs)]
                        for i in range(self.bs)]

            # end_pha = self.smooth(org_pha, win=win)
            end_pha = np.asarray(end_pha)
            bp.setiterblocks(end_pha)

        output = bp.getresult()

        if not iscomplex:
            output = np.angle(output)

        return output

    @staticmethod
    def smooth(array, win):
        return filters.uniform_filter(array.real, win) + 1j * filters.uniform_filter(array.imag, win)

@pyrat.docstringfrom(LFFE_boxcar)
def lffe_boxcar(*args, **kwargs):
    return LFFE_boxcar(*args, **kwargs).run(*args, **kwargs)


class LFFE_Gauss(pyrat.FilterWorker):
    """
    (Local Fringe Frequency Estimator) Gaussian filter with local fringe removal

    :author: Joel Amao
    """
    gui = {'menu': 'InSAR|Phase noise filter', 'entry': 'Local Fringe Gaussian'}
    para = [
        {'var': 'bs', 'value': 32, 'type': 'int', 'text': 'Window size'},
        {'var': 'prc', 'value': 0.5, 'type': 'float', 'text': 'Patch overlap (%)'},
        {'var': 'gDev', 'value': 5.0, 'type': 'float', 'range': [3, 999], 'text': 'Gaussian size'},
        {'var': 'debug', 'value': False, 'type': 'bool', 'text': 'for debugging'},
        {'var': 'thc', 'value': 0.3, 'type': 'float', 'text': 'coherence threshold'},
        {'var': 'gwfun', 'value': False, 'type': 'bool', 'text': 'Use gaussian weights'},
        {'var': 'skipdef', 'value': False, 'type': 'bool', 'text': 'Skips defringing (debug option)'}
    ]

    def __init__(self, *args, **kwargs):
        super(LFFE_Gauss, self).__init__(*args, **kwargs)
        self.name = "LFFE GAUSS PHASE FILTER "
        self.blockprocess = True
        self.blockoverlap = self.bs
        self.scaling_hint = 'phase'
        self.bsx = self.bsy = self.bs
        self.text = 'None'

    def filter(self, array, *args, **kwargs):

        def get_res(block, pha, bs, thc, skipdef):
            # Gets residual (flat) phase and autocorrelation peak
            corr = np.abs(block)
            size = corr.shape[0]
            amax = np.unravel_index(corr.argmax(), corr.shape)
            row = np.fft.fftfreq(size)[amax[0]]
            col = np.fft.fftfreq(size)[amax[1]]
            res = np.zeros_like(pha)

            # b defines the 3-dB value
            coh = np.max(corr)/(bs**2)
            if coh < thc:
                # print("Coherence in patch " + str(coh))
                row = 0
                col = 0

            if skipdef:
                res = pha

            else:
                x = np.linspace(0, self.bs - 1, self.bs) * row
                y = np.linspace(0, self.bs - 1, self.bs) * col
                xv, yv = np.meshgrid(y, x)
                res = pha * np.exp(-1j*2*np.pi*(xv + yv))

            return res, row, col, np.max(corr)

        iscomplex = True
        if not np.iscomplexobj(array):
            # Getting interferogram from phase array
            iscomplex = False
            array = np.exp(1j * array)

        def gaussian_wfun(length):
            if self.bs == 32:
                val = ss.windows.gaussian(length, 5)
            else:
                val = ss.windows.cosine(length)
            return val

        if self.gwfun:
            bp = Blocxy(array, (self.bsx, self.bsy), (self.bsx // int((1 // self.prc)), self.bsy // int((1 // self.prc))),
                    wfunc = gaussian_wfun)
        else:
            bp = Blocxy(array, (self.bsx, self.bsy), (self.bsx // int((1 // self.prc)), self.bsy // int((1 // self.prc))))

        for block in bp.getiterblocks():
            org_pha = block.copy()
            # Interferogram spectrum
            block = np.exp(1j * np.angle(block))
            block = np.fft.fft2(block)

            # Gets residual interferogram
            res_pha, freqx, freqy, peak_val = get_res(block, org_pha, self.bs, self.thc, self.skipdef)

            # Phase smoothing
            fil_pha = self.gauss(res_pha, dev=self.gDev)

            # Adds fringes back
            def add1(arr):
                return np.asarray([[arr[i, j] * np.exp(1j * 2 * np.pi * (freqx * i + freqy * j)) for j in range(self.bs)]
                            for i in range(self.bs)])

            def add2(arr):
                x = np.linspace(0, self.bs - 1, self.bs) * freqx
                y = np.linspace(0, self.bs - 1, self.bs) * freqy
                xv, yv = np.meshgrid(y, x)
                return fil_pha * np.exp(1j*2*np.pi*(xv + yv))

            if self.skipdef:
                end_pha = fil_pha
                # end_pha = add2(fil_pha)
            else:
                end_pha = add2(fil_pha)
            bp.setiterblocks(end_pha)

        output = bp.getresult()

        if not iscomplex:
            output = np.angle(output)

        return output

    @staticmethod
    def gauss(array, dev):
        gLim = 3.0
        return filters.gaussian_filter(array.real, dev, mode='nearest',truncate=gLim) \
               + 1j * filters.gaussian_filter(array.imag, dev, mode='nearest',truncate=gLim)

@pyrat.docstringfrom(LFFE_Gauss)
def lffe_gauss(*args, **kwargs):
    return LFFE_Gauss(*args, **kwargs).run(*args, **kwargs)

class Baran(pyrat.FilterWorker):
    """
    Baran phase noise filter

    further information:
    I.Baran et. al.: A Modification to the Goldstein Radar Interferogram Filer
    IEEE Trans. Geos. and Rem. Sens., Vol. 41, No. 9, pp.2114-2118, 2003


    :author: Andreas Reigber
    """
    gui = {'menu': 'InSAR|Phase noise filter', 'entry': 'Baran'}
    para = [
        {'var': 'bs', 'value': 32, 'type': 'int', 'text': 'Window size'}
    ]

    def __init__(self, *args, **kwargs):
        super(Baran, self).__init__(*args, **kwargs)
        self.name = "BARAN FILTER"
        self.blockprocess = True
        self.blockoverlap = self.bs
        self.scaling_hint = 'phase'

    def filter(self, array, *args, **kwargs):
        array = np.exp(1j * array)
        win = (7, 7)
        coh = filters.uniform_filter(
            np.abs(filters.uniform_filter(array.real, win) + 1j * filters.uniform_filter(array.imag, win)), win)

        bp = Blocxy(array, (self.bs, self.bs), (self.bs // 2, self.bs // 2))
        bc = Blocxy(coh, (self.bs, self.bs), (self.bs // 2, self.bs // 2))
        for block, cblock in zip(bp.getiterblocks(), bc.getiterblocks()):
            block = np.fft.fft2(block)
            block *= abs(block) ** (1.0 - np.mean(cblock[self.bs // 4:self.bs // 4 * 3, self.bs // 4:self.bs // 4 * 3]))
            block = np.fft.ifft2(block)
            bp.setiterblocks(block)
        output = bp.getresult()
        output = np.angle(output)
        return output

@pyrat.docstringfrom(Baran)
def baran(*args, **kwargs):
    return Baran(*args, **kwargs).run(*args, **kwargs)


class QEF(pyrat.FilterWorker):
    """
    QEF phase noise filter. Warning: Memory intensive!

    O. Obsikov et. al.: "Phase noise filtering based on the quasi-ergodic quantum fluctuation theorem",
    MCTV Journal of Modern Quantum Octics, Vol. 17, pp. 4352-4259, 1998.

    :author: Nemjak Nadow
    """
    gui = {'menu': 'InSAR|Phase noise filter', 'entry': 'QEF'}
    para = [
        {'var': 'strength', 'value': 0.5, 'type': 'float', 'range': [0.0, 1.0], 'text': 'Filter strength'},
        {'var': 'win', 'value': 5, 'type': 'int', 'text': 'Window size'}
    ]

    def __init__(self, *args, **kwargs):
        super(QEF, self).__init__(*args, **kwargs)
        self.name = "QEF FILTER"
        self.blockprocess = False
        self.scaling_hint = 'phase'

    def filter(self, array, *args, **kwargs):

        res = np.zeros_like(array)
        res = self.fluc_sub(array, res, 21, 3, 11)
        res = self.fluc_sub(array, res, 9, 3, 11)
        res = self.fluc_sub(array, res, 3, 3, 11)

        res = filters.uniform_filter(res, 11)
        qef = np.exp(1j*(array-res))
        feq = self.median(qef, self.win)
        qef = feq*self.strength + qef*(1-self.strength)
        return np.angle(qef*np.exp(1j*res))

    def fluc_sub(self, pha, upha, a1, alpha1, a2):
        cum = np.exp(1j*(pha-upha))
        if a1 > 1:
            cum = self.smooth(cum, a1)
        if alpha1 > 1:
            cum = np.angle(self.median(cum, alpha1))
        else:
            cum = np.angle(cum)
        ucum = self.ergo_m(cum)
        upha += ucum
        if a2 > 1:
            upha = filters.uniform_filter(upha, a2)
        return upha

    @staticmethod
    def ergo_m(phase):
        ys, xs = phase.shape
        uw = np.angle(np.exp(1j * (np.roll(phase, -1, axis=1) - phase)))
        uw[:, -1] = np.angle(np.exp(1j * (phase[:, -1] - phase[:, -2])))
        x2 = np.angle(np.exp(1j * (phase - np.roll(phase, 1, axis=1))))
        x2[:, 0] = np.angle(np.exp(1j * (phase[:, 1] - phase[:, 0])))
        uw -= x2
        x2 = np.angle(np.exp(1j * (np.roll(phase, -1, axis=0) - phase)))
        x2[-1, :] = np.angle(np.exp(1j * (phase[-1, :] - phase[-2, :])))
        uw += x2
        x2 = np.angle(np.exp(1j * (phase - np.roll(phase, 1, axis=0))))
        x2[0, :] = np.angle(np.exp(1j * (phase[1, :] - phase[0, :])))
        uw -= x2

        for x in range(xs):
            uw[:, x] = np.fft.fft(np.append(uw[:, x], uw[-2:0:-1, x]))[0:ys].real
        for y in range(ys):
            uw[y, :] = np.fft.fft(np.append(uw[y, :], uw[y, -2:0:-1]))[0:xs].real

        wv = [np.cos(np.linspace(0, np.pi, ys)), np.cos(np.linspace(0, np.pi, xs))]
        x2 = reduce(np.add.outer, wv)
        x2[0, 0] = 0.00001

        uw /= (2.0 * (x2 - 2.0))
        uw[0, 0] = 0.0
        for x in range(xs):
            uw[:, x] = np.fft.ifft(np.append(uw[:, x], uw[-2:0:-1, x]))[0:ys].real
        for y in range(ys):
            uw[y, :] = np.fft.ifft(np.append(uw[y, :], uw[y, -2:0:-1]))[0:xs].real
        return uw

    @staticmethod
    def smooth(array, win):
        return filters.uniform_filter(array.real, win) + 1j * filters.uniform_filter(array.imag, win)

    @staticmethod
    def median(array, win):
        return filters.median_filter(array.real, win) + 1j * filters.median_filter(array.imag, win)

@pyrat.docstringfrom(QEF)
def qef(*args, **kwargs):
    return QEF(*args, **kwargs).run(*args, **kwargs)
