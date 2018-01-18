import numpy
from scipy.ndimage import filters

import pyrat
import scipy.fftpack as fftpack
import numpy as np
import logging
import pdb

from pyrat.filter import Unweight


class FFT(pyrat.FilterWorker):
    """
    Fourier transform in range, azimuth or both directions. If multidimensional data are used
    as input, the Fourier transform is applied on each channel separately. The output is centred}
    around the image centre.

    :author: Andreas Reigber
    """

    gui = {'menu': 'SAR|Spectral tools', 'entry': 'Fourier transform'}
    para = [
        {'var': 'axis', 'value': '2D', 'type': 'list', 'range': ['range', 'azimuth', '2D'], 'text': 'FFT axis'},
        {'var': 'method', 'value': 'FFT', 'type': 'list', 'range': ['FFT', 'IFFT'], 'text': 'method'}
    ]

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
        if self.method == 'FFT':
            fft_transform = np.fft.fft
        else:
            fft_transform = np.fft.ifft
        if self.axis == 'range' or self.axis == '2D':
            if self.method == 'IFFT':
                array = np.roll(array, -array.shape[-1] // 2, axis=-1)
            array = fft_transform(array, axis=-1)
            if self.method == 'FFT':
                array = np.roll(array, array.shape[-1] // 2, axis=-1)
        if self.axis == 'azimuth' or self.axis == '2D':
            if self.method == 'IFFT':
                array = np.roll(array, -array.shape[-2] // 2, axis=-2)
            array = fft_transform(array, axis=-2)
            if self.method == 'FFT':
                array = np.roll(array, array.shape[-2] // 2, axis=-2)
        return array


@pyrat.docstringfrom(FFT)
def fft(*args, **kwargs):
    return FFT(*args, **kwargs).run(**kwargs)


class Weight(pyrat.Worker):
    """
    Apply spectral weighting to control sidelobe levels. Image bandwidth is automatically estimated using the cutoff threshold.
    If needed, the image spectrum can be equalised or centred prior to weighting.

    :author: Andreas Reigber
    """
    gui = {'menu': 'SAR|Sidelobe control', 'entry': 'Spectral Weighting'}

    para = [
        {'var': 'func', 'value': 'Hamming', 'type': 'list', 'range': ['Hamming', 'Lanczos'],
         'text': 'Weighting function'},
        {'var': 'alpha', 'value': 0.54, 'type': 'float', 'range': [0.0, 1.0], 'text': 'Weight parameter'},
        {'var': 'axis', 'value': 'both', 'type': 'list', 'range': ['range', 'azimuth', 'both'], 'text': 'Axis'},
        {'var': 'fix', 'value': False, 'type': 'bool', 'text': 'Spectral equalisation'},
        {'var': 'center', 'value': False, 'type': 'bool', 'text': 'Center spectrum'},
        {'var': 'cutoff', 'value': 0.05, 'type': 'float', 'text': 'Cutoff threshold'}
    ]

    def run(self, *args, **kwargs):
        layer = pyrat.data.active

        # STEP1: Estimate spectrum
        self.vblock = False
        rgspec = self.layer_accumulate(self.estimate_spectrum, axis='range', combine=self.combine_spectrum,
                                       silent=False)
        self.vblock = True
        azspec = self.layer_accumulate(self.estimate_spectrum, axis='azimuth', combine=self.combine_spectrum,
                                       silent=False)

        # STEP2: Adjust spectra
        rgcorr, rgcent = self.spec_correction(rgspec, alpha=self.alpha, fix=self.fix, cutoff=self.cutoff,
                                              func=(self.func, self.alpha))
        azcorr, azcent = self.spec_correction(azspec, alpha=self.alpha, fix=self.fix, cutoff=self.cutoff,
                                              func=(self.func, self.alpha))

        if self.center is False:
            azcent, rgcent = None, None

        # STEP3: Weight / Unweight
        if self.axis == 'range' or self.axis == 'both':
            self.vblock = False
            outlayer1 = self.layer_process(self.unweight_spectrum, axis='range', correction=(azcorr, rgcorr),
                                           center=(azcent, rgcent), silent=False)
        if self.axis == 'azimuth' or self.axis == 'both':
            self.vblock = True
            outlayer2 = self.layer_process(self.unweight_spectrum, axis='azimuth', correction=(azcorr, rgcorr),
                                           center=(azcent, rgcent),
                                           layer=outlayer1, silent=False)

        pyrat.delete(outlayer1)
        pyrat.activate(outlayer2)
        return outlayer2

    @staticmethod
    def estimate_spectrum(data, axis='range', **kwargs):
        if axis == 'range':
            spec = np.mean(np.abs(np.fft.fft(data, axis=-1)), axis=-2)
        if axis == 'azimuth':
            spec = np.mean(np.abs(np.fft.fft(data, axis=-2)), axis=-1)
        return spec

    @staticmethod
    def combine_spectrum(values, **kwargs):
        spec = np.zeros_like(values[0])
        for val in values:
            spec += val
        return spec

    @staticmethod
    def spec_correction(array, alpha=1.00, fix=True, cutoff=0.05, func=(None, None)):
        if array.ndim == 1:
            array = array[np.newaxis, ...]
        corr = np.empty_like(array)
        cent = [None] * len(array)
        for k, arr in enumerate(array):
            spec = arr / np.mean(arr)
            spec = filters.uniform_filter(spec, len(array) // 16, mode='wrap')
            offset = np.argmax(spec)
            spec = np.roll(spec, -offset)
            up_bound = 0
            while up_bound < len(spec) and spec[up_bound] > cutoff * spec[0]:
                up_bound += 1
            lo_bound = -1
            while lo_bound > -len(spec) and spec[lo_bound] > cutoff * spec[0]:
                lo_bound -= 1

            if func[0] == "Hamming" or func[0] is None:
                w = func[1] - (1.0 - func[1]) * np.cos(
                    2 * np.pi * np.arange(up_bound - lo_bound) / (up_bound - lo_bound - 1))
            elif func[0] == "Lanczos":
                w = np.sinc(2 * np.arange(up_bound - lo_bound) / (up_bound - lo_bound - 1) - 1)

            if fix is True:
                corr[k, ...] = 1.0 / spec
            else:
                corr[k, ...] = 1.0

            corr[k, up_bound - 1:lo_bound + 1] = 0.0

            corr[k, lo_bound:] *= w[0:-lo_bound]
            corr[k, :up_bound] *= w[-up_bound:]
            corr[k, ...] = np.roll(corr[k, ...], offset)
            corr[k, ...] /= np.mean(corr[k, ...])

            if offset > len(spec) // 2:
                offset -= len(spec)
            cent[k] = (up_bound + lo_bound) // 2 + offset
        return np.squeeze(corr), cent

    @staticmethod
    def unweight_spectrum(data, axis='range', correction=(None, None), center=(None, None), block=(0, 0, 0, 0),
                          **kwargs):
        rgcorr = correction[1]
        azcorr = correction[0]
        # stop()
        if axis == 'range':
            data = np.fft.fft(data, axis=-1)
            data *= rgcorr[..., np.newaxis, :]
            if center[1] is not None:
                data = np.roll(data, -center[1][0], axis=-1)
            data = np.fft.ifft(data, axis=-1)
        if axis == 'azimuth':
            data = np.fft.fft(data, axis=-2)
            data *= azcorr[..., np.newaxis]
            if center[0] is not None:
                data = np.roll(data, -center[0][0], axis=-2)
            data = np.fft.ifft(data, axis=-2)
        return data


@pyrat.docstringfrom(Weight)
def weight(*args, **kwargs):
    return Weight(*args, **kwargs).run(*args, **kwargs)


class Unweight(Weight):
    """
    Spectral unweighting - equalise (and eventually center) the image spectrum. The image bandwidth is
    automatically estimated using the cutoff threshold.

    :author: Andreas Reigber
    """
    gui = {'menu': 'SAR|Sidelobe control', 'entry': 'Unweight spectrum'}

    para = [
        {'var': 'axis', 'value': 'both', 'type': 'list', 'range': ['range', 'azimuth', 'both'], 'text': 'Axis'},
        {'var': 'center', 'value': False, 'type': 'bool', 'text': 'Center spectrum'},
        {'var': 'cutoff', 'value': 0.05, 'type': 'float', 'text': 'Cutoff threshold'}
    ]

    def __init__(self, *args, **kwargs):
        super(Unweight, self).__init__(*args, **kwargs)
        self.name = 'UNWEIGHT'
        self.func = "Hamming"
        self.alpha = 1.0
        self.fix = True


@pyrat.docstringfrom(Unweight)
def unweight(*args, **kwargs):
    return Unweight(*args, **kwargs).run(*args, **kwargs)


class CDA(pyrat.Worker):
    """
    Complex dual apodization filter - a non-linear operation to suppress sidelobes while preserving maximum resolution.
    Can only be used on comlex SLC images and image stacks (PolSAR vector images).

    :author: Andreas Reigber
    """
    gui = {'menu': 'SAR|Sidelobe control', 'entry': 'CDA sidelobe suppression'}

    para = [
        {'var': 'cutoff', 'value': 0.1, 'type': 'float', 'text': 'Cutoff threshold'}
    ]

    def run(self, *args, **kwargs):
        lay0 = pyrat.data.active

        # STEP1: Estimate spectra
        lay1 = pyrat.filter.unweight(layer=lay0, center=True, fix=True, cutoff=self.cutoff)
        lay2 = pyrat.filter.weight(layer=lay0, center=False, fix=False, func='Hamming', alpha=0.5, cutoff=self.cutoff)
        lay3 = self.layer_process(self.cda, layer=[lay1, lay2])

        pyrat.delete(lay1)
        pyrat.delete(lay2)
        pyrat.activate(lay3)

    @staticmethod
    def cda(array, **kwargs):

        c1 = array[0].imag
        c2 = array[1].imag
        im = np.zeros_like(c1)
        aux = (c1 > 0) & (c2 > 0)
        foo = np.minimum(c1, c2)
        im[aux] = foo[aux]
        aux = (c1 < 0) & (c2 < 0)
        foo = np.maximum(c1, c2)
        im[aux] = foo[aux]

        c1 = array[0].real
        c2 = array[1].real
        re = np.zeros_like(c1)
        aux = (c1 > 0) & (c2 > 0)
        foo = np.minimum(c1, c2)
        re[aux] = foo[aux]
        aux = (c1 < 0) & (c2 < 0)
        foo = np.maximum(c1, c2)
        re[aux] = foo[aux]

        return re + 1j * im


@pyrat.docstringfrom(CDA)
def cda(*args, **kwargs):
    return CDA(*args, **kwargs).run(*args, **kwargs)
