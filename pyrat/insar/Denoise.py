import numpy as np
from scipy.ndimage import filters
from functools import reduce

import pyrat
from pyrat.lib.ste import Blocxy
from pyrat.tools import linspace_nd


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
        {'var': 'bs', 'value': 32, 'type': 'int', 'text': 'Window size'}
    ]

    def __init__(self, *args, **kwargs):
        super(Goldstein, self).__init__(*args, **kwargs)
        self.name = "GOLDSTEIN FILTER"
        self.blockprocess = True
        self.blockoverlap = self.bs
        self.scaling_hint = 'phase'

    def filter(self, array, *args, **kwargs):
        array = np.exp(1j * array)
        bp = Blocxy(array, (self.bs, self.bs), (self.bs // 2, self.bs // 2))
        for block in bp.getiterblocks():
            block = np.fft.fft2(block)
            block *= abs(block) ** self.exp
            block = np.fft.ifft2(block)
            bp.setiterblocks(block)
        output = bp.getresult()
        output = np.angle(output)
        return output


@pyrat.docstringfrom(Goldstein)
def goldstein(*args, **kwargs):
    return Goldstein(*args, **kwargs).run(*args, **kwargs)


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
