import pyrat
import numpy as np


class Entalpani(pyrat.FilterWorker):
    """
    Simple Boxcar filter...
    """
    gui = {'menu': 'PolSAR|Parameters', 'entry': 'Entropy / Alpha / Anisotropy'}

    def __init__(self, *args, **kwargs):
        super(Entalpani, self).__init__(*args, **kwargs)    
        self.name = "H/a/A"
        self.blockprocess = True
                
    def filter(self, array, *args, **kwargs):
        
        zdim = array[1].shape[0]
        sew = np.sum(array[0], axis=0)
        pi  = (array[0] + (sew==0)) / (sew + (sew==0))
        
        np.seterr(divide='ignore',invalid='ignore')
        entropy   = np.sum(-pi*np.log(pi)/np.log(zdim), axis=0)
        alphamax  = np.arccos(np.abs(array[1][0,0,...]))
        alphamean = np.sum(np.arccos(np.abs(array[1][0,...]))*pi, axis=0)
        sew = array[0][1,...]+array[0][2,...]
        anisotropy = (array[0][1,...]-array[0][2,...])/(sew+(sew==0))
        np.seterr(divide='warn',invalid='warn')
        return entropy, alphamax, alphamean, anisotropy


def entalpani(*args, **kwargs):
    return Entalpani(*args, **kwargs).run(**kwargs)
