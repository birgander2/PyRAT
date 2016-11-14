import pyrat
import numpy as np
from scipy.ndimage import filters


class Entalpani(pyrat.FilterWorker):
    """
    Estimation of polarimetric entropy, alpha angle (max & mean) and anisotropy.
    Expects a layer with an eigenvalue decomposition as imput.
    """
    gui = {'menu': 'PolSAR|Parameters', 'entry': 'Entropy / Alpha / Anisotropy'}

    def __init__(self, *args, **kwargs):
        super(Entalpani, self).__init__(*args, **kwargs)
        self.name = "H/a/A"
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs["meta"]
        zdim = array[1].shape[0]
        sew = np.sum(array[0], axis=0)
        pi = (array[0] + (sew == 0)) / (sew + (sew == 0))

        np.seterr(divide='ignore', invalid='ignore')
        entropy = np.sum(-pi * np.log(pi) / np.log(zdim), axis=0)
        alphamax = np.arccos(np.abs(array[1][0, 0, ...]))
        alphamean = np.sum(np.arccos(np.abs(array[1][0, ...])) * pi, axis=0)
        sew = array[0][1, ...] + array[0][2, ...]
        anisotropy = (array[0][1, ...] - array[0][2, ...]) / (sew + (sew == 0))
        np.seterr(divide='warn', invalid='warn')
        meta = meta[0]
        meta['CH_name'] = ['Entropy', 'AlphaMax', 'AlphaMean', 'Anisotropy']
        return np.stack([entropy, alphamax, alphamean, anisotropy])


@pyrat.docstringfrom(Entalpani)
def entalpani(*args, **kwargs):
    return Entalpani(*args, **kwargs).run(*args, **kwargs)


class OrientationAngle(pyrat.FilterWorker):
    """
    Estimation of the polarimetric orientation angle.
    """
    gui = {'menu': 'PolSAR|Parameters', 'entry': 'Orientation angle'}
    para = [
        {'var': 'win', 'value': [7, 7], 'type': 'int', 'range': [3, 999], 'text': 'Window size',
         'subtext': ['range', 'azimuth']}
    ]

    def __init__(self, *args, **kwargs):
        super(OrientationAngle, self).__init__(*args, **kwargs)
        self.name = "ORIENTATION ANGLE"
        self.blockprocess = True
        self.blockoverlap = self.win[0] // 2 + 1

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        pol = meta['CH_pol']

        idx_hh = pol.index('HH')
        idx_vv = pol.index('VV')
        idx_xx = pol.index('XX')

        Srr = (array[idx_hh, ...] - array[idx_vv, ...] + 1j * 2.0 * array[idx_xx, ...]) / 2.0
        Sll = (array[idx_vv, ...] - array[idx_hh, ...] + 1j * 2.0 * array[idx_xx, ...]) / 2.0
        a1 = Srr * np.conj(Sll)
        a1 = filters.uniform_filter(a1.real, self.win) + 1j * filters.uniform_filter(a1.imag, self.win)
        return np.angle(a1) / 4


@pyrat.docstringfrom(OrientationAngle)
def orientationangle(*args, **kwargs):
    return OrientationAngle(*args, **kwargs).run(*args, **kwargs)
