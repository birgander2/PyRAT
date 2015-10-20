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
        {'var': 'method', 'value': 'average HV & VH', 'type': 'list', 'range': ['average HV & VH', 'take only HV', 'take only VH'], 'text': 'method :'}
    ]

    def __init__(self, *args, **kwargs):
        super(CalibXsym, self).__init__(*args, **kwargs)
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        pol = meta['CH_pol']

        if array.ndim == 3:                     # PolSAR Vector
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
                oarray[2, ...] = (array[idx_hv, ...] + array[idx_vh, ...])/np.sqrt(2)
            elif self.method == 'take only HV':
                oarray[2, ...] = array[idx_hv, ...]*np.sqrt(2)
            elif self.method == 'take only VH':
                oarray[2, ...] = array[idx_vh, ...]*np.sqrt(2)
            meta['CH_pol'] = ['HH', 'VV', 'XX']
            return oarray
        elif array.ndim == 4:                  # PolSAR Matrix
            idx_hv = pol.index('HVHV*') % 4
            idx_vh = pol.index('VHVH*') % 4
            if self.method == 'average HV & VH':
                # THE FOLLOWING LINE IS BUGGY; needs fixing...
                #array[idx_hv, ...] = (array[idx_hv, ...] + array[idx_vh, ...]) / 2.0
                raise NotImplementedError('TODO!') 
            elif self.method == 'take only HV':
                pass
            elif self.method == 'take only VH':
                 array[idx_hv, ...] =  array[idx_vh, ...]*2
            oarray = np.delete(np.delete(array, idx_vh, axis=0), idx_vh, axis=1)
            return oarray

  
def calibxsym(*args, **kwargs):
    return CalibXsym(*args, **kwargs).run(**kwargs)


class Beta2Gamma(pyrat.FilterWorker):
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

        r0 = meta['rd'] * meta['c0']/2.0 + meta['c0']/meta['rsf']/2.0*np.arange(meta['nrg'])
        h0 = meta['h0']

        for k in range(arr.shape[-2]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                theta = np.arccos((h0 - dem[k, :]) / r0)
            theta[~np.isfinite(theta)] = 0
            if self.type is True:
                arr[..., k, :] *= np.tan(theta)                                          # intensity
            else:
                arr[..., k, :] *= np.sqrt(np.tan(theta))                                 # amplitude

        return arr


def beta2gamma(*args, **kwargs):
    return Beta2Gamma(*args, **kwargs).run(**kwargs)



class NESZFilter(pyrat.FilterWorker):
    """
    Noise floor (NESZ) removal
    """
    gui = {'menu': 'PolSAR|Calibration', 'entry': 'Noise Floor (NESZ) Removal'}
    para = [
        {'var': 'gdev', 'value': 7.0, 'type': 'float', 'range': [1.0, 1e3], 'text': 'Filter Sigma'}
    ]

    def __init__(self, *args, **kwargs):
        super(NESZFilter, self).__init__(*args, **kwargs)
        self.blockprocess = True
        self.blockoverlap = int(3*self.gdev+1.5)
        

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta']
        pol = meta['CH_pol']
        
        win = [self.gdev,self.gdev]
        flt_r = lambda d: filters.gaussian_filter(d, win)
        flt_i = lambda d: flt_r((d*np.conj(d)).real)
        flt_c = lambda d: flt_r(d.real) + 1j*flt_r(d.imag)

        if array.ndim == 3:                     # PolSAR Vector
            idx = [pol.index(p) for p in ['HV','VH']]

            coh = np.abs(flt_c(array[idx[0],:,:]*np.conj(array[idx[1],:,:])))
            coh /= np.sqrt(flt_i(array[idx[0],:,:])*flt_i(array[idx[1],:,:]))

            return array * np.sqrt(coh)[np.newaxis,:,:]
        elif array.ndim == 4:                  # PolSAR Matrix
            idx = [pol.index(p) for p in ['HVHV*','VHVH*','HVVH*']]

            coh = np.abs(flt_c(array[idx[2],:,:]))
            coh /= np.sqrt(flt_r(array[idx[0],:,:].real)*flt_r(array[idx[1],:,:]))
            
            return array * coh[np.newaxis,:,:]

  
def neszfilter(*args, **kwargs):
    return NESZFilter(*args, **kwargs).run(**kwargs)
