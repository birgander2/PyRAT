import PyRat
import numpy as np

import ipdb; stop = ipdb.set_trace

class Cov2Coh(PyRat.FilterWorker):
    """
    Covariance to coherency matrix conversion...
    """

    def __init__(self, *args, **kwargs):
        super(Cov2Coh, self).__init__(*args, **kwargs)    
        self.name = "COVARIANCE -> COHERENCY MATRIX"
        self.allowed_ndim  = [4]
        self.require_para  = ['CH_pol']  
        #self.blockprocess  = True
        #self.nthreads      = 1
        
    def filter(self, array, *args, **kwargs):
        meta  = kwargs['meta']
        pol   = list(meta['CH_pol'])
        if array.shape[0] == 3:
            idx_hh = pol.index('HHHH*') % 3
            idx_vv = pol.index('VVVV*') % 3
            idx_xx = pol.index('XXXX*') % 3
            order = [idx_hh, idx_vv, idx_xx]
            d  = np.array([[1.,1.,0.],[1.,-1.,0.],[0.,0.,np.sqrt(2.)]]) / np.sqrt(2.)
            d = d[order]
            di = np.conj(d).transpose()           
        elif array.shape[0] == 4:
            idx_hh = pol.index('HHHH*') % 4
            idx_vv = pol.index('VVVV*') % 4
            idx_hv = pol.index('HVHV*') % 4
            idx_vh = pol.index('VHVH*') % 4
            d  = np.array([[1.,1.,0.,0.],[1.,-1.,0.,0.],[0.,0.,1.,1.],[0.,0.,-1j,1j]]) / np.sqrt(2)
            di = np.conj(d).transpose()
 
        stop()
        out = d * array * di
        return out
