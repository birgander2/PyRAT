import PyRat
import scipy as sp, numpy as np
import pdb, time

class Eigen(PyRat.FilterWorker):
    """
    Eigendecomposition
    """

    def __init__(self, *args, **kwargs):
        super(Eigen, self).__init__(*args, **kwargs)    
        self.name = "EIGEN DECOMPOSITION"
        self.allowed_ndim  = [4]
        self.blockprocess = True
                
    def filter(self, array, *args, **kwargs):
        attrs = kwargs['meta']

        vdim, zdim = array.shape[0:2]
        ydim, xdim = array.shape[2:4]
        newshape = (vdim+1, zdim, ydim, xdim)
        
        array  = np.rollaxis(np.rollaxis(array,0,start=4),0,start=4)
        ew, ev = np.linalg.eigh(array)
        ew[ew <= 0] = 0.0
        sidx   = np.indices(ew.shape)
        idx    = np.argsort(-ew,axis=2)
        ew     = ew[sidx[0],sidx[1],idx]
        
        sidx   = np.indices(ev.shape)
        for k in range(zdim):
            sidx[3,:,:,k,:] = idx
        ev = ev[sidx[0],sidx[1],sidx[2],sidx[3]]

        del attrs['CH_pol']
        return np.rollaxis(ew,2), np.rollaxis(np.rollaxis(ev,3),3)
  
