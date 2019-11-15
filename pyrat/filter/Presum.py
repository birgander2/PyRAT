import pyrat
from pyrat.filter.Filter import Median
from pyrat.filter.tools import rebin
import PIL
import numpy as np


class Presum(pyrat.Worker):
    """
    Presumming of data arrays. Improved version working 'inplace', i.e. without reading the entire data set.

    :author: Andreas Reigber
    """
    gui = {'menu': 'Tools', 'entry': 'Presumming'}
    para = [
        {'var': 'subx', 'value': 4, 'type': 'int', 'range': [0, 999], 'text': 'Presumming range'},
        {'var': 'suby', 'value': 4, 'type': 'int', 'range': [0, 999], 'text': 'Presumming azimuth'},
        {'var': 'decimate', 'value': False, 'type': 'bool', 'text': 'Skip averaging'}
    ]

    def __init__(self, *args, **kwargs):
        super(Presum, self).__init__(*args, **kwargs)
        self.name = "PRESUMMING"
        self.blockprocess = False
        self.nthreads = 1

    def run(self, *args, **kwargs):
        meta = pyrat.data.getAnnotation(layer=self.layer)
        li = pyrat.query(layer=self.layer)
        odim = li.shape
        nry = odim[-2]

        odim[-2] //= self.suby
        odim[-1] //= self.subx

        outlayer = pyrat.data.addLayer(dtype=li.dtype, shape=odim)
        blockdim = odim.copy()
        blockdim[-2] = 1

        P = pyrat.tools.ProgressBar('  ' + self.name, odim[-2])
        P.update(0)
        for k in range(odim[-2]):
            arr = pyrat.getdata(block=(k*self.suby, (k+1)*self.suby, 0, odim[-1] * self.subx), layer=self.layer)
            if self.decimate is True:
                arr = arr[..., ::self.suby, ::self.subx]
            else:
                arr = rebin(arr, tuple(blockdim))
            pyrat.data.setData(arr, block=(k, k+1, 0, 0), layer=outlayer)
        P.update(k + 1)
        del P

        pyrat.activate(outlayer)

        if "geo_ps_east" in meta and "geo_ps_north" in meta:
            meta['geo_ps_east'] = meta["geo_ps_east"] * self.subx
            meta['geo_ps_north'] = meta["geo_ps_north"] * self.suby
        pyrat.data.setAnnotation(meta, layer=outlayer)
        return outlayer


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
        {'var': 'subx', 'value': 4, 'type': 'int', 'range': [0, 999], 'text': 'Lanczos range'},
        {'var': 'suby', 'value': 4, 'type': 'int', 'range': [0, 999], 'text': 'Lanczos azimuth'},
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
            arr = array[k, ...]  # workaround for PIL bug
            while True:  # (large array crash)
                try:
                    arr = PIL.Image.fromarray(arr.astype(np.float32), mode='F')
                except OverflowError:
                    shp = arr.shape
                    arr = rebin(arr[0:shp[-2] // 2 * 2, 0:shp[-1] // 2 * 2], (shp[-2] // 2, shp[-1] // 2))
                else:
                    break

            # arr = PIL.Image.fromarray(array[k, ...].astype(np.float32), mode='F')
            arr = arr.resize(odim[-2:][::-1], resample=method)
            oarr[k, ...] = np.asarray(arr)

        return np.squeeze(oarr)


@pyrat.docstringfrom(Lanczos)
def lanczos(*args, **kwargs):
    return Lanczos(*args, **kwargs).run(**kwargs)


@pyrat.docstringfrom(Median)
def medianf(*args, **kwargs):
    return Median(*args, **kwargs).run(**kwargs)