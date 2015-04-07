import pyrat
import numpy as np


class Lex2Pauli(pyrat.FilterWorker):
    """
    Lexicographic to Pauli conversion...
    """
    gui = {'menu': 'PolSAR|Decompositions', 'entry': 'Pauli decomposition'}

    def __init__(self, *args, **kwargs):
        super(Lex2Pauli, self).__init__(*args, **kwargs)
        self.name = "LEXICOGRAPHIC TO PAULI CONVERSION"
        self.allowed_ndim = [3]
        self.require_para = ['CH_pol']
        self.blockprocess = True
        # self.nthreads      = 1

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        pol = meta['CH_pol']
        oarray = np.empty_like(array)
        self.x = 1.7
        if array.shape[0] == 3:
            idx_hh = pol.index('HH')
            idx_vv = pol.index('VV')
            idx_xx = pol.index('XX')
            oarray[0, ...] = array[idx_hh, ...] + array[idx_vv, ...]
            oarray[1, ...] = array[idx_hh, ...] - array[idx_vv, ...]
            oarray[2, ...] = np.sqrt(2) * array[idx_xx, ...]
            oarray /= np.sqrt(2)
            meta['CH_pol'] = ['P1', 'P2', 'P3']
        elif array.shape[0] == 4:
            idx_hh = pol.index('HH')
            idx_vv = pol.index('VV')
            idx_hv = pol.index('HV')
            idx_vh = pol.index('VH')
            oarray[0, ...] = array[idx_hh, ...] + array[idx_vv, ...]
            oarray[1, ...] = array[idx_hh, ...] - array[idx_vv, ...]
            oarray[2, ...] = array[idx_hv, ...] + array[idx_vh, ...]
            oarray[3, ...] = array[idx_vh, ...] - array[idx_hv, ...]
            oarray /= np.sqrt(2)
            meta['CH_pol'] = ['P1', 'P2', 'P3', 'P4']
        return oarray


def lex2pauli(*args, **kwargs):
    return Lex2Pauli(*args, **kwargs).run(**kwargs)


class Eigen(pyrat.FilterWorker):
    """
    Eigendecomposition
    """
    gui = {'menu': 'PolSAR|Decompositions', 'entry': 'Eigenvector / Eigenvalue'}

    def __init__(self, *args, **kwargs):
        super(Eigen, self).__init__(*args, **kwargs)
        self.name = "EIGEN DECOMPOSITION"
        self.allowed_ndim = [4]
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        attrs = kwargs['meta']

        vdim, zdim = array.shape[0:2]
        ydim, xdim = array.shape[2:4]
        newshape = (vdim + 1, zdim, ydim, xdim)

        array = np.rollaxis(np.rollaxis(array, 0, start=4), 0, start=4)
        ew, ev = np.linalg.eigh(array)
        ew[ew <= 0] = 0.0
        sidx = np.indices(ew.shape)
        idx = np.argsort(-ew, axis=2)
        ew = ew[sidx[0], sidx[1], idx]

        sidx = np.indices(ev.shape)
        for k in range(zdim):
            sidx[3, :, :, k, :] = idx
        ev = ev[sidx[0], sidx[1], sidx[2], sidx[3]]

        del attrs['CH_pol']
        return np.rollaxis(ew, 2), np.rollaxis(np.rollaxis(ev, 3), 3)


def eigen(*args, **kwargs):
    return Eigen(*args, **kwargs).run(**kwargs)
