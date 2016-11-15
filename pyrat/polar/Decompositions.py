import numpy as np
import logging
import pyrat
from .tools import C_to_T

# try:
#     # from .Cy_Descompositions import cy_yamaguchi4
#
#     class Yamaguchi4(pyrat.FilterWorker):
#         """
#
#         :author: Andreas Reigber
#         """
#         gui = {'menu': 'PolSAR|Decompositions', 'entry': 'Yamaguchi4'}
#
#         def __init__(self, *args, **kwargs):
#             pyrat.FilterWorker.__init__(self, *args, **kwargs)
#             self.name = "YAMAGUCHI4"
#             self.allowed_ndim = [4]
#             self.require_para = ['CH_pol']
#             self.blockprocess = False
#
#         def filter(self, array, *args, **kwargs):
#             meta = kwargs["meta"]
#             pol = meta['CH_pol']
#
#             array = C_to_T(array, pol)                           # covariance to coherency
#             shp = array.shape
#             array = np.reshape(array, (shp[0] * shp[1], shp[2], shp[3]))
#             idx_p1 = pol.index('P1P1*')
#             idx_p2 = pol.index('P2P2*')
#             idx_p3 = pol.index('P3P3*')
#             if len(pol) == 9:                                       # T3 Matrix
#                 cal_factor = 2.0
#                 npol = 3
#             else:                                                   # T4 Matrix
#                 cal_factor = 1.0
#                 npol = 4
#             span = np.abs(np.trace(np.reshape(array, (npol, npol, array.shape[-2], array.shape[-1]))))
#
#             np.seterr(divide='ignore', invalid='ignore')
#
#
#
#
#             array = cy_yamaguchi4(array)
#
#             pd = np.clip(fd ), 0.0, span)                             # clip to valid range
#             ps = np.clip(fs ), 0.0, span)                             # clip to valid range
#             pv = np.clip(fv, 0.0, span)                               # clip to valid range
#             ph = np.clip(fh, 0.0, span)                               # clip to valid range
#
#             np.seterr(divide='warn', invalid='warn')                  # remove nans and zeros
#             array = np.stack([pd, ps, pv, ph])
#             array[~np.isfinite(array)] = 0.0
#             array[array < 0.0] = 0.0
#
#             del meta['CH_pol']
#             meta['CH_name'] = ['even', 'odd', 'volume', 'helix']
#             return array
#
#
#     @pyrat.docstringfrom(Yamaguchi4)
#     def yamaguchi4(*args, **kwargs):
#         return Yamaguchi4(*args, **kwargs).run(*args, **kwargs)
# except ImportError:
#     logging.info("Yamaguchi4 cython module not found. (run build process?)")


class Yamaguchi3(pyrat.FilterWorker):
    """
    Yamaguchi 3-component decomposition into surface, double bounce and volume component.
    Expects as input a C3 or C4 polarimetric covariance matrix. Returns the three estimated power
    components in the order even-odd-volume.

    :author: Andreas Reigber
    """
    gui = {'menu': 'PolSAR|Decompositions', 'entry': 'Yamaguchi3'}

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.name = "YAMAGUCHI3"
        self.allowed_ndim = [4]
        self.require_para = ['CH_pol']
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs["meta"]
        pol = meta['CH_pol']

        np.seterr(divide='ignore', invalid='ignore')

        shp = array.shape
        array = np.reshape(array, (shp[0] * shp[1], shp[2], shp[3]))

        idx_hh = pol.index('HHHH*')
        idx_vv = pol.index('VVVV*')
        idx_hhvv = pol.index('HHVV*')
        if len(pol) == 9:                                       # C3 Matrix
            idx_xx = pol.index('XXXX*')
            cal_factor = 2.0
            npol = 3
        else:                                                   # C4 Matrix
            idx_xx = pol.index('HVHV*')
            cal_factor = 1.0
            npol = 4

        span = np.abs(np.trace(np.reshape(array, (npol, npol, array.shape[-2], array.shape[-1]))))

        hhhh = np.abs(array[idx_hh, ...])
        vvvv = np.abs(array[idx_vv, ...])
        hvhv = np.abs(array[idx_xx, ...]) / cal_factor
        hhvvr = np.array(np.abs(array[idx_hhvv, ...]).real)
        hhvvi = np.array(np.abs(array[idx_hhvv, ...]).imag)

        pd, ps, pv = self.decomp_yamaguchi3(hhhh, vvvv, hvhv, hhvvr, hhvvi, span)

        np.seterr(divide='warn', invalid='warn')                # remove nans and zeros
        array = np.stack([pd, ps, pv])
        array[~np.isfinite(array)] = 0.0
        array[array < 0.0] = 0.0

        del meta['CH_pol']
        meta['CH_name'] = ['even', 'odd', 'volume']
        return array

    @staticmethod
    def decomp_yamaguchi3(hhhh, vvvv, hvhv, hhvvr, hhvvi, span):

        fv = np.zeros_like(hhhh)
        fd = np.zeros_like(hhhh)
        fs = np.zeros_like(hhhh)
        alp_r = np.zeros_like(hhhh)
        alp_i = np.zeros_like(hhhh)
        bet_r = np.zeros_like(hhhh)
        bet_i = np.zeros_like(hhhh)

        ratio = 10*np.log10(vvvv / hhhh)                       # preconditioning

        idx_m2 = ratio < -2.0
        fv[idx_m2] = 15.0 * hvhv[idx_m2] / 4.0
        hhhh[idx_m2] -= 8.0 * fv[idx_m2] / 15.0
        vvvv[idx_m2] -= 3.0 * fv[idx_m2] / 15.0
        hhvvr[idx_m2] -= 2.0 * fv[idx_m2] / 15.0

        idx_p2 = ratio > 2.0
        fv[idx_p2] = 15.0 * hvhv[idx_p2] / 4.0
        hhhh[idx_p2] -= 3.0 * fv[idx_p2] / 15.0
        vvvv[idx_p2] -= 8.0 * fv[idx_p2] / 15.0
        hhvvr[idx_p2] -= 2.0 * fv[idx_p2] / 15.0

        idx_mm = (ratio >= -2.0) & (ratio <= 2.0)
        fv[idx_mm] = 8.0 * hvhv[idx_mm] / 2.0
        hhhh[idx_mm] -= 3.0 * fv[idx_mm] / 8.0
        vvvv[idx_mm] -= 3.0 * fv[idx_mm] / 8.0
        hhvvr[idx_mm] -= 1.0 * fv[idx_mm] / 8.0

        idx_volerr = (hhhh <= 0.0) | (vvvv <= 0.0)                # volume > total
        fv[idx_volerr] = span[idx_volerr] # - hvhv[idx_volerr]

        foo = hhvvr**2 + hhvvi**2                                # non-realisable Shhvv term
        idx = (foo > hhhh*vvvv) & ~idx_volerr
        hhvvr[idx] *= np.sqrt(hhhh[idx]*vvvv[idx]/foo[idx])
        hhvvi[idx] *= np.sqrt(hhhh[idx]*vvvv[idx]/foo[idx])

        idx = (hhvvr >= 0.0) & ~idx_volerr
        alp_r[idx] = -1.0                                            # single dominant
        alp_i[idx] = 0.0
        fd[idx] = (hhhh[idx] * vvvv[idx] - hhvvr[idx]**2 - hhvvi[idx]**2) / (hhhh[idx] + vvvv[idx] + 2 * hhvvr[idx])
        fs[idx] = vvvv[idx] - fd[idx]
        bet_r[idx] = (fd[idx] + hhvvr[idx]) / fs[idx]
        bet_i[idx] = hhvvi[idx] / fs[idx]

        idx = (hhvvr < 0.0) & ~idx_volerr
        bet_r[idx] = 1.0                                                  # double dominant
        bet_i[idx] = 0.0
        fs[idx] = (hhhh[idx] * vvvv[idx] - hhvvr[idx]**2 - hhvvi[idx]**2) / (hhhh[idx] + vvvv[idx] - 2 * hhvvr[idx])
        fd[idx] = vvvv[idx] - fs[idx]
        alp_r[idx] = (hhvvr[idx] - fs[idx]) / fd[idx]
        alp_i[idx] = hhvvi[idx] / fd[idx]

        pd = np.clip(fd * (1 + alp_r**2 + alp_i**2), 0.0, span)                             # clip to valid range
        ps = np.clip(fs * (1 + bet_r**2 + bet_i**2), 0.0, span)                             # clip to valid range
        pv = np.clip(fv, 0.0, span)                                                         # clip to valid range

        return pd, ps, pv


@pyrat.docstringfrom(Yamaguchi3)
def yamaguchi3(*args, **kwargs):
    return Yamaguchi3(*args, **kwargs).run(**kwargs)


class FreemanDurden(pyrat.FilterWorker):
    """
    Freeman Durden decomposition into surface, double bounce and volume component. Expects
    as input a C3 or C4 polarimetric covariance matrix. Returns the three estimated power
    components in the order even-odd-volume.

    A. Freeman, S.L. Durden: "A three-component scattering model for polarimetric SAR data"
    Transactions on Geoscience and Remote Sensing, Vol. 36. No. 3, pp. 963-973, 1998

    :author: Andreas Reigber
    :tested: for C4 input (todo: C3)
    """
    gui = {'menu': 'PolSAR|Decompositions', 'entry': 'Freeman-Durden'}

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.name = "FreemanDurden"
        self.allowed_ndim = [4]
        self.require_para = ['CH_pol']
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs["meta"]
        pol = meta['CH_pol']

        np.seterr(divide='ignore', invalid='ignore')

        shp = array.shape
        array = np.reshape(array, (shp[0] * shp[1], shp[2], shp[3]))

        idx_hh = pol.index('HHHH*')
        idx_vv = pol.index('VVVV*')
        idx_hhvv = pol.index('HHVV*')
        if len(pol) == 9:                                       # C3 Matrix
            idx_xx = pol.index('XXXX*')
            cal_factor = 2.0
            npol = 3
        else:                                                   # C4 Matrix
            idx_xx = pol.index('HVHV*')
            cal_factor = 1.0
            npol = 4

        span = np.abs(np.trace(np.reshape(array, (npol, npol, array.shape[-2], array.shape[-1]))))

        fv = np.abs(array[idx_xx, ...]) * 3 / cal_factor        # volume component
        pv = 8 * fv / 3

        pd = np.zeros_like(pv)
        ps = np.zeros_like(pv)
        shh = np.abs(array[idx_hh, ...]) - fv
        svv = np.abs(array[idx_vv, ...]) - fv
        shv = array[idx_hhvv, ...] - fv / 3

        idx_volerr = (shh <= 0.0) | (svv <= 0.0)                # volume > total
        pv[idx_volerr] = span[idx_volerr]

        foo = np.abs(shv)**2                                    # data conditioning for
        idx = (foo > shh*svv) & ~idx_volerr                     # non-realisable ShhShv term
        shv[idx] *= np.sqrt(shh[idx]*svv[idx]/foo[idx])

        idx = (shv.real >= 0.0) & ~idx_volerr
        alpha = -1.0                                             # double dominant
        c1 = shh[idx]
        c2 = svv[idx]
        c3r = shv[idx].real
        c3i = shv[idx].imag
        fd = (c1 * c2 - c3r**2 - c3i**2) / (c1 + c2 + 2 * c3r)
        fs = c2 - fd
        beta = np.sqrt((fd + c3r)**2 + c3i**2) / fs
        ps[idx] = fs * (1 + beta**2)
        pd[idx] = fd * (1 + alpha**2)

        idx = (shv.real < 0.0) & ~idx_volerr
        beta = 1.0                                              # surface dominant
        c1 = shh[idx]
        c2 = svv[idx]
        c3r = shv[idx].real
        c3i = shv[idx].imag
        fs = (c1 * c2 - c3r ** 2 - c3i ** 2) / (c1 + c2 - 2 * c3r)
        fd = c2 - fs
        alpha = np.sqrt((fs - c3r) ** 2 + c3i ** 2) / fd
        ps[idx] = fs * (1 + beta ** 2)
        pd[idx] = fd * (1 + alpha ** 2)

        pd = np.clip(pd, 0.0, span)                             # clip to valid range
        ps = np.clip(ps, 0.0, span)
        pv = np.clip(pv, 0.0, span)

        np.seterr(divide='warn', invalid='warn')                # remove nan's
        array = np.stack([pd, ps, pv])
        array[~np.isfinite(array)] = 0.0
        array[array < 0.0] = 0.0

        del meta['CH_pol']
        meta['CH_name'] = ['even', 'odd', 'volume']
        return array


@pyrat.docstringfrom(FreemanDurden)
def freemandurden(*args, **kwargs):
    return FreemanDurden(*args, **kwargs).run(**kwargs)


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

        # todo: replace with static method from PolSAR worker
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
