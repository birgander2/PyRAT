import numpy as np
import itertools
from scipy.ndimage import filters


def smooth(array, box, phase=False):
    """
    Imitates IDL's smooth function. Can also (correctly) smooth interferometric phases with the phase=True keyword.
    """
    if np.iscomplexobj(array):
        return filters.uniform_filter(array.real, box) + 1j * filters.uniform_filter(array.imag, box)
    elif phase is True:
        return np.angle(smooth(np.exp(1j * array), box))
    else:
        return filters.uniform_filter(array.real, box)


def rebin(arr, *shape, **kwargs):
    """
    Imitates IDL's rebin function. Allows also phase rebining.

    :author: Andreas Reigber
    """

    phase = False  # combination of *args and fixed keywords
    if 'phase' in kwargs:
        phase = kwargs['phase']  # works only in in python 3!

    if len(shape) == 1:  # allows to pass shape as list/array or normal arguments
        if type(shape[0]) == int:
            shape = [shape[0]]
        else:
            shape = list(shape[0])

    oarr = arr.copy()
    oshap = oarr.shape
    for d in range(arr.ndim):
        n1 = shape[d]
        n2 = oshap[d]
        if n1 < n2:
            s = list(oarr.shape)
            s.insert(d + 1, n2 // n1)
            s[d] = n1
            if phase is True:
                oarr = np.angle(np.exp(1j * oarr.reshape(s)).mean(d + 1))
            else:
                oarr = oarr.reshape(s).mean(d + 1).astype(oarr.dtype)
        elif n1 > n2:
            oarr = oarr.repeat(n1 // n2, axis=d)
        else:
            pass
    return oarr


def polyfit2d(x, y, z, order=3):

    ncols = (order + 1) ** 2
    g = np.zeros((x.size, ncols))
    ij = itertools.product(range(order + 1), range(order + 1))
    for k, (i, j) in enumerate(ij):
        g[:, k] = x ** i * y ** j
    m, _, _, _ = np.linalg.lstsq(g, z)
    return m


def polyval2d(x, y, m):
    # x = x * 1.0
    #y = y * 1.0
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order + 1), range(order + 1))
    z = np.zeros_like(x, dtype='float64')
    for a, (i, j) in zip(m, ij):
        z += a * x ** i * y ** j
    return z


def coreg(arr1, arr2, sub=False):
    """
    Returns the coregistration offset between two complex arrays

    :author: Andreas Reigber
    :param arr1: The master image as (2D np.ndarray)
    :type arr1: complex
    :param arr2: The slave image as (2D np.ndarray)
    :type arr2: complex
    :param sub=True: Calculate subpixel offset.
    :type arr2: boolean
    :returns: tuple containing slave image offset in y and x
    """
    s = arr1.shape
    out1 = np.fft.rfft2(abs(arr1)).astype('complex64')
    out2 = np.fft.rfft2(abs(arr2)).astype('complex64')
    out1 *= np.conj(out2)
    out2 = np.fft.irfft2(out1).astype('float32')
    pos = out2.argmax()
    xoff = pos % s[1]
    yoff = pos / s[1]
    if xoff >= s[1] // 2:
        xoff = xoff - s[1]
    if yoff >= s[0] // 2:
        yoff = yoff - s[0]
    if sub is True:  # explicit polynominal expansion of the peak
        ysub = 0.5 * (out2[yoff - 1, xoff] - out2[yoff + 1, xoff]) / (out2[yoff - 1, xoff] +
                                                                      out2[yoff + 1, xoff] - 2 * out2[yoff, xoff])
        xsub = 0.5 * (out2[yoff, xoff - 1] - out2[yoff, xoff + 1]) / (out2[yoff, xoff - 1] +
                                                                      out2[yoff, xoff + 1] - 2 * out2[yoff, xoff])
        yoff += np.round(ysub, 2)
        xoff += np.round(xsub, 2)
    return yoff.astype('float32'), xoff.astype('float32')
