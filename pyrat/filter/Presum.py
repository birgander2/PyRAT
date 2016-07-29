from __future__ import print_function
import pyrat
from pyrat.filter.tools import rebin
import PIL
import numpy as np

class Presum(pyrat.FilterWorker):
    """
    Presumming of data arrays.

    :author: Andreas Reigber
    """
    gui = {'menu': 'Tools', 'entry': 'Presumming'}
    para = [
        {'var': 'subx', 'value': 4, 'type': 'int',  'range': [0, 999], 'text': 'Presumming range'},
        {'var': 'suby', 'value': 4, 'type': 'int',  'range': [0, 999], 'text': 'Presumming azimuth'}
        ]

    def __init__(self, *args, **kwargs):
        super(Presum, self).__init__(*args, **kwargs)
        self.name = "PRESUMMING"
        self.blockprocess = False

    def filter(self, array, *args, **kwargs):

        odim = list(array.shape)
        odim[-2] //= self.suby
        odim[-1] //= self.subx

        array = rebin(array[..., 0:odim[-2] * self.suby, 0:odim[-1] * self.subx], odim)

        return array


@pyrat.docstringfrom(Presum)
def presum(*args, **kwargs):
    return Presum(*args, **kwargs).run(**kwargs)


class Lanczos(pyrat.FilterWorker):
    """
    Lanczos-based presumming of data arrays.
    Lanczos is a sinc-type resampling method and should provide better quality than simple rebinning.
    WARNING: Only for float32 arrays & not for matrix images!!!!

    :author: Andreas Reigber
    """

    para = [
        {'var': 'subx', 'value': 4, 'type': 'int',  'range': [0, 999], 'text': 'Lanczos range'},
        {'var': 'suby', 'value': 4, 'type': 'int',  'range': [0, 999], 'text': 'Lanczos azimuth'},
        {'var': 'method', 'value': True, 'type': 'bool'}
        ]

    def __init__(self, *args, **kwargs):
        super(Lanczos, self).__init__(*args, **kwargs)
        self.name = "LANCZOS"
        self.blockprocess = False

    def filter(self, array, *args, **kwargs):
        if self.method is True:
            method = PIL.Image.LANCZOS
        else:
            method = PIL.Image.NEAREST

        odim = list(array.shape)
        odim[-2] //= self.suby
        odim[-1] //= self.subx

        oarr = np.empty(odim, dtype='float32')
        if oarr.ndim == 3:
            nchannel = odim[0]
        else:
            nchannel = 1
            array = array[None, ...]
            oarr = oarr[None, ...]

        for k in range(nchannel):
            arr = array[k, ...]                                                # workaround for PIL bug
            while True:                                                        # (large array crash)
                try:
                    arr = PIL.Image.fromarray(arr.astype(np.float32), mode='F')
                except OverflowError:
                    shp = arr.shape
                    arr = rebin(arr[0:shp[-2]//2*2, 0:shp[-1]//2*2], (shp[-2]//2, shp[-1]//2))
                else:
                    break

            # arr = PIL.Image.fromarray(array[k, ...].astype(np.float32), mode='F')
            arr = arr.resize(odim[-2:][::-1], resample=method)
            oarr[k, ...] = np.asarray(arr)

        return np.squeeze(oarr)


@pyrat.docstringfrom(Lanczos)
def lanczos(*args, **kwargs):
    return Lanczos(*args, **kwargs).run(**kwargs)

