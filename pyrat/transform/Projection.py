import numpy as np
import pyrat
from pyrat.lib.ste import interpol_cubic
from pyrat.lib.ste import interpol_lin1d



class Slant2Ground(pyrat.FilterWorker):
    """
    Slant to ground range projection

    :author: Andreas Reigber
    """
    gui = {'menu': 'SAR|Geometry', 'entry': 'Slant->Ground range'}

    para = [
        {'var': 'minangle', 'value': 20.0, 'type': 'float', 'range': [1.0, 90.0], 'text': 'Minimum near range angle'}
    ]

    def __init__(self, *args, **kwargs):
        super(Slant2Ground, self).__init__(*args, **kwargs)
        self.name = "SLANT2GROUND"
        self.blockprocess = True
        self.require_para = ["h0", "terrain", "rd", "nrg", "rsf"]

    def filter(self, array, *args, **kwargs):
        meta = kwargs["meta"]

        if array.dtype == 'uint8':
            convert = array.dtype
            array = array.astype(np.float32)
        else:
            convert = False

        shp = array.shape
        if array.ndim > 2:
            nchannels = int(np.prod(shp[0:-2]))
            array = array.reshape((nchannels, shp[-2], shp[-1]))
        else:
            nchannels = 1
            array = array[np.newaxis, ...]

        height = meta['h0'] - meta['terrain']
        sr = meta['rd'] * meta['c0'] / 2.0 + np.arange(meta['nrg']) * meta['c0'] / meta['rsf'] / 2.0

        theta = np.arccos(height / sr)
        ang_mask = (theta >= np.radians(self.minangle))
        theta = theta[ang_mask]

        gr = np.sqrt(sr ** 2 - height ** 2)
        gr_near = np.min(gr[ang_mask])
        gr_far = np.max(gr[ang_mask])
        ps_gr = meta['ps_rg'] / np.max(np.sin(theta))

        nout = int((gr_far - gr_near) / ps_gr)
        shp_out = list(array.shape)
        shp_out[-1] = nout
        out = np.empty(shp_out, dtype=array.dtype)

        tin = np.arange(meta['nrg'])
        gr_out = gr_near + np.arange(nout) * ps_gr
        tout = np.float32(interpol_lin1d(gr, tin, gr_out))

        for ch in range(nchannels):
            for x in range(array.shape[-2]):
                out[ch, x, :] = interpol_cubic(array[ch, x, :], tout, threads=1)

        meta['ps_rg'] = ps_gr
        if convert is not False:
            out[out < 0.5] = 0.0
            out[out >= 0.5] = 1.0
            out = out.astype(np.uint8)

        out[~np.isfinite(out)] = 0

        return out


@pyrat.docstringfrom(Slant2Ground)
def slant2ground(*args, **kwargs):
    return Slant2Ground(*args, **kwargs).run(*args, **kwargs)
