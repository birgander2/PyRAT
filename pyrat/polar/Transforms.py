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
        self.blockprocess = False

    def filter(self, array, *args, **kwargs):
        attrs = kwargs['meta']
        chpol = attrs['CH_pol']
        ch_idx = [i for i, x in enumerate(chpol) if x == self.pol.upper()]
        attrs['CH_pol'] =[chpol[ch_idx[k]] for k in range(np.size(ch_idx))]
        return np.squeeze(array[ch_idx,...])


@pyrat.docstringfrom(PolExtract)
def polextract(*args, **kwargs):
    return PolExtract(*args, **kwargs).run(*args, **kwargs)


class PolTrace(pyrat.FilterWorker):
    """
    Extract the trace of a polarimetric matrix

    :author: Marc Jaeger
    """
    gui = {'menu': 'PolSAR|Transform', 'entry': 'Matrix Trace'}
    para = [
        {'var': 'average', 'value': True, 'type': 'bool', 'text': 'average trace by number of channels'}    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "Matrix Trace"
        self.allowed_ndim = [4]
        self.blockprocess = False

    def filter(self, array, *args, **kwargs):
        attrs = kwargs['meta']
        attrs['CH_pol'] = attrs['CH_pol'][0][:2]
        return np.einsum('ii...', array).real / (array.shape[0] if self.average else 1)


@pyrat.docstringfrom(PolTrace)
def poltrace(*args, **kwargs):
    return PolTrace(*args, **kwargs).run(*args, **kwargs)



class Quad2Compact(pyrat.FilterWorker):
    """
    Transforms quadpol data into compact pol dual-channel data

    :author: Andreas Reigber
    """
    gui = {'menu': 'PolSAR|Transform', 'entry': 'Quad -> Compact'}
    para = [
        {'var': 'mode', 'value': 'CTLR', 'type': 'list', 'range': ['CTLR', 'PI4', 'DCP'], 'text': 'Compact mode'}]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "CTLR"
        self.allowed_ndim = [3]
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        attrs = kwargs['meta']
        pol = attrs['CH_pol']
        idx_hh = pol.index('HH')
        idx_vv = pol.index('VV')
        if len(pol) == 4:
            idx_hv = pol.index('HV')
            idx_vh = pol.index('VH')
        elif len(pol) == 3:
            idx_hv = pol.index('XX')
            idx_vh = pol.index('XX')

        shp = list(array.shape)
        shp[0] = 2
        oarray = np.empty(shp, dtype=array.dtype)
        if self.mode == "CTLR":
            oarray[0, ...] = 1.0 / np.sqrt(2) * (array[idx_hh, ...] - 1j * array[idx_hv, ...])
            oarray[1, ...] = 1.0 / np.sqrt(2) * (array[idx_vh, ...] - 1j * array[idx_vv, ...])
            attrs['CH_pol'] = ["CTLR_1", "CTLR_2"]
        if self.mode == "PI4":
            oarray[0, ...] = 1.0 / np.sqrt(2) * (array[idx_hh, ...] + array[idx_hv, ...])
            oarray[1, ...] = 1.0 / np.sqrt(2) * (array[idx_vh, ...] + array[idx_vv, ...])
            attrs['CH_pol'] = ["PI4_1", "PI4_2"]
        if self.mode == "DCP":
            oarray[0, ...] = 0.5 * (array[idx_hh, ...] - array[idx_vv, ...] + 1j * array[idx_hv, ...] + 1j * array[idx_vh, ...])
            oarray[1, ...] = 0.5 * (1j * array[idx_hh, ...] + 1j * array[idx_vv, ...] + array[idx_vh, ...] - array[idx_hv, ...])
            attrs['CH_pol'] = ["DCP_1", "DCP_2"]
        return oarray


@pyrat.docstringfrom(Quad2Compact)
def quad2compact(*args, **kwargs):
    return Quad2Compact(*args, **kwargs).run(*args, **kwargs)

