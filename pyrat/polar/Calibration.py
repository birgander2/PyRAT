import pyrat
import numpy as np


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
