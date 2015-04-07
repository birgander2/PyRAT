import pyrat
import numpy as np


class Vec2Mat(pyrat.FilterWorker):
    """
    Vector to Matrix transform...

    :author: Andreas Reigber
    """
    gui = {'menu': 'PolSAR|Transform', 'entry': 'Vector -> Matrix'}

    def __init__(self, *args, **kwargs):
        super(Vec2Mat, self).__init__(*args, **kwargs)    
        self.name = "VECTOR -> MATRIX"       
        self.allowed_ndim = [3]
        self.blockprocess = True
        
    def filter(self, array, *args, **kwargs):
        out = array[np.newaxis, :, :, :] * array[:, np.newaxis, :, :].conjugate()
        attrs = kwargs['meta']
        if 'CH_pol' in attrs:
            attrs['CH_pol'] = [p1+p2+'*' for p1 in attrs['CH_pol'] for p2 in attrs['CH_pol']]
        return out


def vec2mat(*args, **kwargs):
    return Vec2Mat(*args, **kwargs).run(**kwargs)
