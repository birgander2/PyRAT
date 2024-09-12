import pyrat
import numpy as np
from scipy.ndimage import filters
import warnings

class CalibXsym(pyrat.FilterWorker):
    """
    Cross-polar symmetrisation
    """
    gui = {'menu': 'PolSAR|Calibration', 'entry': 'Cross-polar symmetrisation'}
    para = [
        {'var': 'method', 'value': 'average HV & VH', 'type': 'list',
         'range': ['average HV & VH', 'take only HV', 'take only VH'], 'text': 'method :'}
    ]

    def __init__(self, *args, **kwargs):
        super(CalibXsym, self).__init__(*args, **kwargs)
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        pol = meta['CH_pol']
        if array.ndim == 3:  # PolSAR Vector
            idx_hh = pol.index('HH')
            idx_vv = pol.index('VV')
            idx_hv = pol.index('HV')
            idx_vh = pol.index('VH')
            shp = list(array.shape)
            shp[0] = 3
            oarray = np.empty(shp, dtype=array.dtype)
            oarray[0, ...] = array[idx_hh, ...]
            oarray[1, ...] = array[idx_vv, ...]
            if self.method == 'average HV & VH':
                oarray[2, ...] = (array[idx_hv, ...] + array[idx_vh, ...]) / np.sqrt(2)
            elif self.method == 'take only HV':
                oarray[2, ...] = array[idx_hv, ...] * np.sqrt(2)
            elif self.method == 'take only VH':
                oarray[2, ...] = array[idx_vh, ...] * np.sqrt(2)
            meta['CH_pol'] = ['HH', 'VV', 'XX']
            return oarray
        elif array.ndim == 4:  # PolSAR Matrix
            idx_hvhv = pol.index('HVHV*')
            idx_vhhv = pol.index('VHHV*')
            idx_hvvh = pol.index('HVVH*')
            idx_vhvh = pol.index('VHVH*')
            idx_hhhv = pol.index('HHHV*')
            idx_hhvh = pol.index('HHVH*')
            idx_vvhv = pol.index('VVHV*')
            idx_vvvh = pol.index('VVVH*')
            idx_hvhh = pol.index('HVHH*')
            idx_vhhh = pol.index('VHHH*')
            idx_hvvv = pol.index('HVVV*')
            idx_vhvv = pol.index('VHVV*')
            shp = array.shape
            array = np.reshape(array, (shp[0] * shp[1], shp[2], shp[3]))
            mask = np.ones(len(array), dtype=bool)
            if self.method == 'average HV & VH':
                arr_hhhv = (array[idx_hhhv, ...] + array[idx_hhvh, ...]) / np.sqrt(2.0)
                arr_vvhv = (array[idx_vvhv, ...] + array[idx_vvvh, ...]) / np.sqrt(2.0)
                arr_xxxx = (array[idx_hvhv, ...] + array[idx_vhvh, ...])
                array[idx_hhhv, ...] = arr_hhhv
                array[idx_hvhh, ...] = np.conj(arr_hhhv)
                array[idx_vvhv, ...] = arr_vvhv
                array[idx_hvvv, ...] = np.conj(arr_vvhv)
                array[idx_hvhv, ...] = arr_xxxx
                mask[[idx_hhvh, idx_vhhh, idx_vvvh, idx_vhvv, idx_hvvh, idx_vhhv, idx_vhvh]] = False
            elif self.method == 'take only HV':
                array = np.delete(array, [idx_hhvh, idx_vvvh, idx_hvvh, idx_vhvh, idx_vhhh, idx_vhvv, idx_vhhv], axis=0)
            elif self.method == 'take only VH':
                array = np.delete(array, [idx_hhhv, idx_vvhv, idx_hvhv, idx_vhhv, idx_hvhh, idx_hvvv, idx_hvvh], axis=0)
            array = array[mask, ...]
            pol = [p for p, b in zip(pol, mask) if b == True]
            for k, p in enumerate(pol):
                if p[0:2] in ['HV', 'VH']:
                    pol[k] = 'XX' + pol[k][2:]
                if p[2:4] in ['HV', 'VH']:
                    pol[k] = pol[k][0:2] + 'XX*'
            array = np.reshape(array, (shp[0] - 1, shp[1] - 1, shp[2], shp[3]))
            meta['CH_pol'] = pol
            return array


@pyrat.docstringfrom(CalibXsym)
def calibxsym(*args, **kwargs):
    return CalibXsym(*args, **kwargs).run(*args, **kwargs)

class Beta2Gamma(pyrat.FilterWorker):
    """
    Beta0 to Gamma0 conversion
    """
    para = [
        {'var': 'type', 'value': False, 'type': 'bool', 'text': 'intensity / covariance'}
    ]

    def __init__(self, *args, **kwargs):
        super(Beta2Gamma, self).__init__(*args, **kwargs)
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta'][0]
        arr = array[0]
        dem = array[1].astype(np.float32)

        r0 = meta['rd'] * meta['c0'] / 2.0 + meta['c0'] / meta['rsf'] / 2.0 * np.arange(meta['nrg'])
        h0 = meta['h0']

        for k in range(arr.shape[-2]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                theta = np.arccos((h0 - dem[k, :]) / r0)
            theta[~np.isfinite(theta)] = 0
            if self.type is True:
                arr[..., k, :] *= np.tan(theta)  # intensity
            else:
                arr[..., k, :] *= np.sqrt(np.tan(theta))  # amplitude

        return arr


@pyrat.docstringfrom(Beta2Gamma)
def beta2gamma(*args, **kwargs):
    return Beta2Gamma(*args, **kwargs).run(*args, **kwargs)

class Beta2Sigma(pyrat.FilterWorker):
    """
    Beta0 to Sigma0 conversion
    """
    para = [
        {'var': 'type', 'value': False, 'type': 'bool', 'text': 'intensity / covariance'}
    ]

    def __init__(self, *args, **kwargs):
        super(Beta2Sigma, self).__init__(*args, **kwargs)
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta'][0]
        arr = array[0]
        dem = array[1].astype(np.float32)

        r0 = meta['rd'] * meta['c0'] / 2.0 + meta['c0'] / meta['rsf'] / 2.0 * np.arange(meta['nrg'])
        h0 = meta['h0']

        for k in range(arr.shape[-2]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                theta = np.arccos((h0 - dem[k, :]) / r0)
            theta[~np.isfinite(theta)] = 0
            if self.type is True:
                arr[..., k, :] *= np.sin(theta)  # intensity
            else:
                arr[..., k, :] *= np.sqrt(np.sin(theta))  # amplitude

        return arr


@pyrat.docstringfrom(Beta2Sigma)
def beta2sigma(*args, **kwargs):
    return Beta2Sigma(*args, **kwargs).run(*args, **kwargs)


class NESZFilter(pyrat.FilterWorker):
    """
    Noise floor (NESZ) removal
    """
    gui = {'menu': 'PolSAR|Calibration', 'entry': 'Noise Floor (NESZ) Removal'}
    para = [
        {'var': 'gdev', 'value': 13.0, 'type': 'float', 'range': [1.0, 1e3], 'text': 'Filter Sigma'},
        {'var': 'noise', 'value': 1.0, 'type': 'float', 'range': [0.0, 2.0], 'text': 'Noise removal factor'}
    ]

    def __init__(self, *args, **kwargs):
        super(NESZFilter, self).__init__(*args, **kwargs)
        self.blockprocess = True
        self.blockoverlap = int(3 * self.gdev + 1.5)

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        pol = meta['CH_pol']

        win = [self.gdev, self.gdev]
        flt_r = lambda d: filters.gaussian_filter(d, win)
        flt_i = lambda d: flt_r((d * np.conj(d)).real)
        flt_c = lambda d: flt_r(d.real) + 1j * flt_r(d.imag)

        if array.ndim == 3:  # PolSAR Vector
            idx = [pol.index(p) for p in ['HV', 'VH']]

            sig = np.abs(flt_c(array[idx[0], ...] * np.conj(array[idx[1], ...])))
            signoise = np.sqrt(flt_i(array[idx[0], ...]) * flt_i(array[idx[1], ...]))
            noise = signoise - sig
            np.seterr(divide='ignore', invalid='ignore')
            for k in range(len(array)):
                intensity = abs(array[k, ...]) ** 2
                array[k, ...] *= np.sqrt((intensity - noise - self.noise * np.sqrt(intensity * noise)) / intensity)
            array[~np.isfinite(array)] = 0.0
            np.seterr(divide='warn', invalid='warn')
            return array
            # return array * np.sqrt(coh)[np.newaxis, :, :]
        elif array.ndim == 4:  # PolSAR Matrix
            idx = [pol.index(p) for p in ['HVHV*', 'VHVH*', 'HVVH*']]

            coh = np.abs(flt_c(array[idx[2], :, :]))
            coh /= np.sqrt(flt_r(array[idx[0], :, :].real) * flt_r(array[idx[1], :, :]))
            print('WARNING: Matrix NESZ removal buggy - Copolar intensities not considered!')
            return array * coh[np.newaxis, :, :]


@pyrat.docstringfrom(NESZFilter)
def neszfilter(*args, **kwargs):
    return NESZFilter(*args, **kwargs).run(**kwargs)