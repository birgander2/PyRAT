import pyrat
import numpy as np


class Lex2Pauli(pyrat.FilterWorker):
    """
    Lexicographic to Pauli basis conversion.
    This routine works both on lexicographic scattering vectors andcovariance matrices as input. Covariance matices
    are converted to coherency matrices. Input matrices can be 3x3 or 4x4. Compact polarimetric modes are not supported.
    """
    gui = {'menu': 'PolSAR|Decompositions', 'entry': 'Pauli decomposition'}

    def __init__(self, *args, **kwargs):
        super(Lex2Pauli, self).__init__(*args, **kwargs)
        self.name = "LEXICOGRAPHIC TO PAULI CONVERSION"
        self.allowed_ndim = [3, 4]
        self.require_para = ['CH_pol']
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        pol = meta['CH_pol']
        oarray = np.empty_like(array)
        if array.ndim == 3:  # vector transform
            if array.shape[0] == 3:
                idx_hh = pol.index('HH')
                idx_vv = pol.index('VV')
                idx_xx = pol.index('XX')
                oarray[0, ...] = (array[idx_hh, ...] + array[idx_vv, ...])
                oarray[1, ...] = (array[idx_hh, ...] - array[idx_vv, ...])
                oarray[2, ...] = array[idx_xx, ...] * np.sqrt(2)
                oarray /= np.sqrt(2)
                meta['CH_pol'] = ['P1', 'P2', 'P3']
            elif array.shape[0] == 4:
                idx_hh = pol.index('HH')
                idx_vv = pol.index('VV')
                idx_hv = pol.index('HV')
                idx_vh = pol.index('VH')
                oarray[0, ...] = (array[idx_hh, ...] + array[idx_vv, ...])
                oarray[1, ...] = (array[idx_hh, ...] - array[idx_vv, ...])
                oarray[2, ...] = (array[idx_hv, ...] + array[idx_vh, ...])
                oarray[3, ...] = (array[idx_hv, ...] - array[idx_vh, ...]) * 1j
                oarray /= np.sqrt(2)
                meta['CH_pol'] = ['P1', 'P2', 'P3', 'P4']
        if array.ndim == 4:  # matrix transform
            if array.shape[0] == 3:
                idx_hh = pol.index('HHHH*') % 3
                idx_vv = pol.index('VVVV*') % 3
                idx_xx = pol.index('XXXX*') % 3
                idx = np.array([idx_hh, idx_vv, idx_xx])
                D = np.array([[1.0, 1.0, 0.0],
                              [1.0, -1.0, 0.0],
                              [0.0, 0.0, np.sqrt(2)]], dtype='f4') / np.sqrt(2.0)
                oarray = array[idx, ...][:, idx, ...]
                oarray = np.rollaxis(np.rollaxis(oarray, 0, start=4), 0, start=4)
                oarray = np.matmul(D, np.matmul(oarray, np.conj(D.T)))
                oarray = np.rollaxis(np.rollaxis(oarray, 3), 3)
                meta['CH_pol'] = ['P1P1*', 'P1P2*', 'P1P3*', 'P2P1*', 'P2P2*', 'P2P3*', 'P3P1*', 'P3P2*', 'P3P3*']
            elif array.shape[0] == 4:
                idx_hh = pol.index('HHHH*') % 4
                idx_vv = pol.index('VVVV*') % 4
                idx_hv = pol.index('HVHV*') % 4
                idx_vh = pol.index('VHVH*') % 4
                idx = np.array([idx_hh, idx_vv, idx_hv, idx_vh])
                D = np.array([[1.0, 1.0, 0.0, 0.0],
                              [1.0, -1.0, 0.0, 0.0],
                              [0.0, 0.0, 1.0, 1.0],
                              [0.0, 0.0, -1j, 1j]], dtype='c8') / np.sqrt(2.0)
                oarray = array[idx, ...][:, idx, ...]
                oarray = np.rollaxis(np.rollaxis(oarray, 0, start=4), 0, start=4)
                oarray = np.matmul(D, np.matmul(oarray, np.conj(D.T)))
                oarray = np.rollaxis(np.rollaxis(oarray, 3), 3)
                meta['CH_pol'] = ['P1P1*', 'P1P2*', 'P1P3*', 'P1P4*', 'P2P1*', 'P2P2*', 'P2P3*', 'P2P4*',
                                  'P3P1*', 'P3P2*', 'P3P3*', 'P3P4*', 'P4P1*', 'P4P2*', 'P4P3*', 'P4P4*']
        return oarray


@pyrat.docstringfrom(Lex2Pauli)
def lex2pauli(*args, **kwargs):
    return Lex2Pauli(*args, **kwargs).run(*args, **kwargs)


class Eigen(pyrat.FilterWorker):
    """
    Eigendecomposition a la Cloude / Pottier.
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


@pyrat.docstringfrom(Eigen)
def eigen(*args, **kwargs):
    return Eigen(*args, **kwargs).run(*args, **kwargs)
