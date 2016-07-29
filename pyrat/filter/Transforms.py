import pyrat
import numpy as np


class Complex2Abs(pyrat.FilterWorker):
    """
    Complex to absolute conversion
    """
    gui = {'menu': 'SAR|Transform', 'entry': 'SLC -> Amplitude'}

    def __init__(self, *args, **kwargs):
        super(Complex2Abs, self).__init__(*args, **kwargs)
        self.name = "Complex -> Abs"
        self.blockprocess = True
        
    def filter(self, array, *args, **kwargs):
        return np.abs(array)


class Int2Amp(pyrat.FilterWorker):
    """
    Intensity to amplitude conversion
    """
    gui = {'menu': 'SAR|Transform', 'entry': 'Intensity -> Amplitude'}

    def __init__(self, *args, **kwargs):
        super(Int2Amp, self).__init__(*args, **kwargs)
        self.name = "INT -> AMP"
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        return np.sqrt(array)


class Amp2Int(pyrat.FilterWorker):
    """
    Amplitude to intensity conversion
    """
    gui = {'menu': 'SAR|Transform', 'entry': 'Amplitude -> Intensity'}

    def __init__(self, *args, **kwargs):
        super(Amp2Int, self).__init__(*args, **kwargs)
        self.name = "AMP -> INT"
        self.blockprocess = True

    def filter(self, array, *args, **kwargs):
        return array**2


@pyrat.docstringfrom(Complex2Abs)
def complex2abs(*args, **kwargs):
    return Complex2Abs(*args, **kwargs).run(**kwargs)


@pyrat.docstringfrom(Amp2Int)
def amp2int(*args, **kwargs):
    return Amp2Int(*args, **kwargs).run(**kwargs)


@pyrat.docstringfrom(Int2Amp)
def int2amp(*args, **kwargs):
    return Int2Amp(*args, **kwargs).run(**kwargs)
