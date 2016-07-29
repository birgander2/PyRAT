from __future__ import print_function
import pyrat
import scipy.fftpack as fftpack
import numpy as np
import logging


class FFT(pyrat.FilterWorker):
    """
    Fourier transform in range, azimuth or both directions. If multidimensional data are used
    as input, the Fourier transform is applied on each channel separately. The output is centred
    around the image centre.

    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Spectral tools', 'entry': 'Fourier transform'}
    para = [{'var': 'axis', 'value': '2D', 'type': 'list', 'range': ['range', 'azimuth', '2D'], 'text': 'FFT axis'}]

    def __init__(self, *args, **kwargs):
        super(FFT, self).__init__(*args, **kwargs)
        self.name = "FOURIER TRANSFORM"
        if self.axis == 'range':
            self.blockprocess = True
        elif self.axis == 'azimuth':
            self.blockprocess = True
            self.vblock = True
        else:
            self.blockprocess = False

    def filter(self, array, *args, **kwargs):

        if self.axis == 'range':
            oarray = np.roll(fftpack.fft(array, axis=-1), array.shape[-1]//2, axis=-1)
        elif self.axis == 'azimuth':
            oarray = np.roll(fftpack.fft(array, axis=-2), array.shape[-2]//2, axis=-2)
        else:
            oarray = np.roll(np.roll(fftpack.fft2(array), array.shape[-1]//2, axis=-1), array.shape[-2]//2, axis=-2)
        return oarray


@pyrat.docstringfrom(FFT)
def fft(*args, **kwargs):
    return FFT(*args, **kwargs).run(**kwargs)
