import PyRat
import numpy as np
import pdb, logging

class Lex2Pauli(PyRat.FilterWorker):
    """
    Lexicographic to Pauli conversion...
    """

    def __init__(self, *args, **kwargs):
        super(Lex2Pauli, self).__init__(*args, **kwargs)    
        self.name = "LEXICOGRAPHIC TO PAULI CONVERSION"
        self.allowed_ndim  = [3]
        self.require_para  = ['CH_pol']  
        self.blockprocess  = True
        #self.nthreads      = 1
        
    def filter(self, array, *args, **kwargs):
        meta  = kwargs['meta']
        pol   = meta['CH_pol']
        oarray = np.empty_like(array)
        self.x = 1.7
        if array.shape[0] == 3:
            idx_hh = pol.index('HH')
            idx_vv = pol.index('VV')
            idx_xx = pol.index('XX')
            oarray[0,...] = array[idx_hh,...] + array[idx_vv,...]
            oarray[1,...] = array[idx_hh,...] - array[idx_vv,...]
            oarray[2,...] = np.sqrt(2)*array[idx_xx,...]
            oarray /= np.sqrt(2)
            meta['CH_pol'] = ['P1','P2','P3']
        elif array.shape[0] == 4:
            idx_hh = pol.index('HH')
            idx_vv = pol.index('VV')
            idx_hv = pol.index('HV')
            idx_vh = pol.index('VH')
            oarray[0,...] = array[idx_hh,...] + array[idx_vv,...]
            oarray[1,...] = array[idx_hh,...] - array[idx_vv,...]
            oarray[2,...] = array[idx_hv,...] + array[idx_vh,...]
            oarray[3,...] = array[idx_vh,...] - array[idx_hv,...]
            oarray /= np.sqrt(2)
            meta['CH_pol'] = ['P1','P2','P3','P4']
        return oarray
