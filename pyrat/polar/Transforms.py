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
        attrs = kwargs['meta']
        if 'CH_pol' in attrs:  # sort a bit the channels if possible
            pol = attrs['CH_pol']
            if len(set(pol) & set(["HH", "VV", "XX"])) == 3:
                idx = [pol.index('HH'), pol.index('VV'), pol.index('XX')]
                array = array[idx, ...]
                attrs['CH_pol'] = [pol[i] for i in idx]
            if len(set(pol) & set(["HH", "VV", "HV", "VH"])) == 4:
                idx = [pol.index('HH'), pol.index('VV'), pol.index('HV'), pol.index('VH')]
                array = array[idx, ...]
                attrs['CH_pol'] = [pol[i] for i in idx]

        out = array[np.newaxis, :, :, :] * array[:, np.newaxis, :, :].conjugate()
        if 'CH_pol' in attrs:
            attrs['CH_pol'] = [p1 + p2 + '*' for p1 in attrs['CH_pol'] for p2 in attrs['CH_pol']]
        return out


@pyrat.docstringfrom(Vec2Mat)
def vec2mat(*args, **kwargs):
    return Vec2Mat(*args, **kwargs).run(*args, **kwargs)


class Rotate(pyrat.FilterWorker):
    """
    Polarimetric basis rotation (pixel by pixel).
    Input are two layers: One with the PolSAR scattering vector in lexicographic basis,
    and a second one containing the rotation angles in radians (possibly estimated
    by pyrat.polar.orientationangle).
    """

    def __init__(self, *args, **kwargs):
        super(Rotate, self).__init__(*args, **kwargs)
        self.name = "POLSAR ROTATION"
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        meta = kwargs['meta'][0]
        pol = meta['CH_pol']

        data = array[0]
        shp = data.shape
        angle = array[1]
        arr = np.empty((shp[1], shp[2], 2, 2), dtype=data.dtype)

        idx_hh = pol.index('HH')
        arr[..., 0, 0] = data[idx_hh, ...]
        idx_vv = pol.index('VV')
        arr[..., 1, 1] = data[idx_vv, ...]
        if 'XX' in pol:
            idx_xx = pol.index('XX')
            arr[..., 0, 1] = data[idx_xx, ...]
            arr[..., 1, 0] = data[idx_xx, ...]
        else:
            idx_hv = pol.index('HV')
            idx_vh = pol.index('VH')
            arr[..., 0, 1] = data[idx_hv, ...]
            arr[..., 1, 0] = data[idx_vh, ...]

        R1 = np.empty((shp[1], shp[2], 2, 2), dtype=angle.dtype)
        R2 = np.empty((shp[1], shp[2], 2, 2), dtype=angle.dtype)
        cosangle = np.cos(angle)
        sinangle = np.sin(angle)
        R1[..., 0, 0] = cosangle
        R1[..., 1, 1] = cosangle
        R1[..., 0, 1] = -sinangle
        R1[..., 1, 0] = sinangle
        R2[..., 0, 0] = cosangle
        R2[..., 1, 1] = cosangle
        R2[..., 0, 1] = sinangle
        R2[..., 1, 0] = -sinangle

        arr = np.matmul(R2, np.matmul(arr, R1))
        out = np.empty_like(data)
        out[idx_hh, ...] = arr[..., 0, 0]
        out[idx_vv, ...] = arr[..., 1, 1]
        if 'XX' in pol:
            out[idx_xx, ...] = arr[..., 1, 0]
        else:
            out[idx_hv, ...] = arr[..., 0, 1]
            out[idx_vh, ...] = arr[..., 1, 0]
        return out


@pyrat.docstringfrom(Rotate)
def rotate(*args, **kwargs):
    return Rotate(*args, **kwargs).run(*args, **kwargs)



class PolExtract(pyrat.FilterWorker):
    """
    Extract a single polarimetric channel

    :author: Marc Jaeger
    """
    gui = {'menu': 'PolSAR|Transform', 'entry': 'Extract Channel'}
    para = [
        {'var': 'pol', 'value': 'HH', 'type': 'list', 'range': ['HH','HV','VH','VV','XX'], 'text': 'Polarisation'}
    ]

    def __init__(self, *args, **kwargs):
        super(PolExtract, self).__init__(*args, **kwargs)
        self.name = "Extract Channel"
        self.allowed_ndim = [3]
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        attrs = kwargs['meta']
        ch_idx = attrs['CH_pol'].index(self.pol)
        attrs['CH_pol'] = self.pol
        return array[ch_idx,...]

