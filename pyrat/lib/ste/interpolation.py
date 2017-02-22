import scipy as sp
from scipy import interpolate
import numpy as np


def interpol_spline1d(x, y, xout):
    """
    Performs a cubic spline interpolation of an irregularily sampled 1D signal

    :author: Andreas Reigber
    :param x: The abscissa values of the signal
    :type x: 1-D ndarray float
    :param y: The ordinate values of the signal
    :type y: 1D ndarray
    :param xout: The abscissa values where the interpolates are desired
    :type xout: 1D ndarray float
    :returns: The interpolated values
    """
    f = interpolate.interp1d(x, y, kind='cubic')
    return f(xout)


def interpol_spline2d(x, y, z, xout, yout):
    """
    Performs a cubic spline interpolation of a 2D matrix/image irregularily sampled on both axes

    :author: Andreas Reigber
    :param x: The x values of the first axis of z
    :type x: 1-D ndarray float
    :param y: The y values of the second axis of z
    :type y: 1-D ndarray float
    :param z: The input matrix
    :type z: 2D ndarray
    :param xout: The values on the first axis where the interpolates are desired
    :type xout: 1D ndarray float
    :param yout: The values on the second axis where the interpolates are desired
    :type yout: 1D ndarray float
    :returns: The interpolated matrix /  image
    """
    if np.iscomplexobj(z):
        f_real = interpolate.interp2d(x, y, z.real, kind='cubic')
        f_imag = interpolate.interp2d(x, y, z.imag, kind='cubic')
        return f_real(xout, yout) + 1j * f_imag(xout, yout)
    else:
        f = interpolate.interp2d(x, y, z, kind='cubic')
        return f(xout, yout)


def interpol_cubic(y, xi, **kwargs):
    """
    1D cubic convolution interpolation on regularly gridded input. Cython version is
    multithreaded and supports a 'threads' keyword to set the number of threads
    (default: # of processors)

    :author: Andreas Reigber
    :param y: The function values
    :type y: 1-D ndarray float
    :param xi: The positions where the interpolates are desired
    :type xi: 1-D ndarray float

    :returns: The interpolated signal
    """
    n = len(y)
    yi = np.empty(xi.shape, y.dtype)
    for i in range(len(xi)):
        if xi[i] < 0.0 or xi[i] > n-1:
            print('WARNING: Bad x input to cubiconv ==> 0 <= x <= len(y)')
            yi[i] = 0.0
        else:
            klo = int(xi[i])
            khi = klo + 1

            if float(klo) == xi[i]:
                yi[i] = y[klo]
            else:
                kmi = klo - 1
                kpl = khi + 1

                if kmi < 0:
                    a = 3 * y[klo] - 3 * y[khi] + y[kpl]
                else:
                    a = y[kmi]

                if kpl > n-1:
                    d = 3 * y[khi] - 3 * y[klo] + y[kmi]
                else:
                    d = y[kpl]

                b = y[klo]
                c = y[khi]
                t = (xi[i] - klo)
                t2 = t * t
                t3 = t2 * t
                c00 = (- t3 + 2 * t2 - t) / 2.0
                c10 = (3 * t3 - 5 * t2 + 2) / 2.0
                c20 = (- 3 * t3 + 4 * t2 + t) / 2.0
                c30 = (t3 - t2) / 2.0
                yi[i] = a * c00 + b * c10 + c * c20 + d * c30
    return yi


def interpol_cubic_irr(x, y, xi, sort=True, **kwargs):
    """
    1D cubic convolution interpolation on irregularly gridded input. Cython version is
    multithreaded and supports a 'threads' keyword to set the number of threads
    (default: # of processors)

    :author: Andreas Reigber
    :param x: The abscissa values
    :type y: 1-D ndarray float
    :param y: The ordinate values
    :type y: 1-D ndarray float, same length as x
    :param xi: The positions where the interpolates are desired
    :type xi: 1-D ndarray float

    :returns: The interpolated signal
    """

    if sort is True:
        sidx = np.argsort(x)
        x = x[sidx]
        y = y[sidx]

    n = len(y)
    yi = np.empty(xi.shape, y.dtype)
    for i in range(len(xi)):
        if xi[i] < x[0] or xi[i] > x[-1]:
            print('WARNING: Bad x input to cubiconv ==> 0 <= x <= len(y)')
            yi[i] = np.nan

        klo = 0
        khi = n-1
        while khi - klo > 1:
            k = int((khi + klo) / 2.0)
            if x[k] > xi[i]:
                khi = k
            else:
                klo = k
        h = x[khi] - x[klo]
        if h == 0.0:
            print('WARNING: Bad x input to cubiconv ==> x values must be distinct')
            yi[i] = np.nan

        kmi = klo - 1
        kpl = khi + 1

        if kmi < 0:
            a = 3 * y[klo] - 3 * y[khi] + y[kpl]
        else:
            a = y[kmi]

        if kpl > n-1:
            d = 3 * y[khi] - 3 * y[klo] + y[kmi]
        else:
            d = y[kpl]

        b = y[klo]
        c = y[khi]
        t = (xi[i] - x[klo]) / h
        t2 = t * t
        t3 = t2 * t
        h2 = h * h
        c00 = (- t3 + 2 * t2 - t) / 2.0
        c10 = (3 * t3 - 5 * t2 + 2) / 2.0
        c20 = (- 3 * t3 + 4 * t2 + t) / 2.0
        c30 = (t3 - t2) / 2.0
        yi[i] = a * c00 + b * c10 + c * c20 + d * c30
    return yi


def interpol_lanczos(y, xi, a=2, **kwargs):
    """
    1D lanczos interpolation on regularly gridded input. Cython version is
    multithreaded and supports a 'threads' keyword to set the number of threads
    (default: # of processors)

    :author: Andreas Reigber
    :param y: The function values
    :type y: 1-D ndarray float
    :param xi: The positions where the interpolates are desired
    :type xi: 1-D ndarray float
    :param a: Lanczos filter size (default=2)
    :type a: int

    :returns: The interpolated signal
    """
    n = len(y)
    yi = np.empty(xi.shape, y.dtype)
    if a < 1:
        a = 1
    a2 = 2*a + 1
    if n < a2:
        print("ERROR: window size too large")
        return
    for i in range(len(xi)):
        yi[i] = 0.0
        if xi[i] < 0.0 or xi[i] > n-1:
            yi[i] = np.nan
            continue
        x1 = int(xi[i]) - a + 1
        if x1 < 0:
            x1 = 0
        if x1 + a2 - 1 > n:
            x1 = n - a2 + 1
        x2 = x1 + a2 - 1
        for x in range(x1, x2):
            if np.abs(xi[i] - x) < a:
                w = np.sinc((xi[i] - x)/a)
            else:
                w = 0
            yi[i] += y[x] * np.sinc((xi[i] - x)) * w
    return yi


def interpol_sinc(y, xi, trunc=0, win='rect', **kwargs):
    """
    1D sinc interpolation on regularly gridded input. No Cython version yet - use
    interpolate_lanczos if performance is needed...

    :author: Andreas Reigber
    :param y: The function values
    :type y: 1-D ndarray float
    :param xi: The positions where the interpolates are desired
    :type xi: 1-D ndarray float
    :param trunc: window length (0=all samples)
    :type trunc: int
    :param win: windowing function ('rect','lanczos','sqcos','hamming','blackman')
    :type win: string

    :returns: The interpolated signal
    """
    n = len(y)
    yi = np.empty(xi.shape, y.dtype)
    if n < trunc:
        print("ERROR: window size too large")
        return
    if trunc == 0:
        trunc = n
    c2 = trunc // 2
    for i in range(len(xi)):
        yi[i] = 0
        if xi[i] < 0.0 or xi[i] > n-1:
            yi[i] = np.nan
            continue
        if trunc == n:
            x1 = 0
            x2 = n
        else:
            x1 = int(xi[i]) - c2 + 1
            if x1 < 0:
                x1 = 0
            if x1 + trunc > n:
                x1 = n - trunc
            x2 = x1 + trunc
        for x in range(x1, x2):
            if np.abs(xi[i] - x) > c2:
                w = 0.0
            else:
                if win == 'blackman':
                    alpha = 0.16
                    a0 = (1.0 - alpha) / 2.0
                    a1 = 0.5
                    a2 = alpha / 2.0
                    w = a0 + a1 * np.cos(2 * np.pi * (xi[i] - x) / trunc) + a2 * np.cos(4 * np.pi * (xi[i] - x) / trunc)
                elif win == 'sqcos':
                    w = np.cos(np.pi * (xi[i] - x) / trunc) * np.cos(np.pi * (xi[i] - x) / trunc)
                elif win == 'hamming':
                    w = 0.54 + 0.46 * np.cos(2 * np.pi * (xi[i] - x) / trunc)
                elif win == 'lanczos':
                    w = np.sinc((xi[i] - x) / c2)
                else:
                    w = 1.0
            yi[i] += y[x] * np.sinc((xi[i] - x)) * w
    return yi


def interpol2D_lanczos(array, yi, xi, a=2, **kwargs):
    """
    2D lanczos interpolation on regularly gridded input. Cython version is
    multithreaded and supports a 'threads' keyword to set the number of threads
    (default: # of processors)

    :author: Andreas Reigber
    :param array: The input array
    :type array: 2-D ndarray float / complex
    :param yi: The y sample positions where the interpolates are desired
    :type yi: 1-D ndarray float
    :param xi: The x sample positions where the interpolates are desired
    :type xi: 1-D ndarray float
    :param a: Lanczos filter size parameter (default=2)
    :type a: int

    :returns: The interpolated signal
    """
    ny = array.shape[0]
    nx = array.shape[1]
    out = np.empty((len(yi), len(xi)), array.dtype)
    if a < 1:
        a = 1
    a2 = 2*a + 1
    if ny < a or nx < a:
        print("ERROR: window size too large")
        return

    for xp in range(len(xi)):
        for yp in range(len(yi)):
            out[yp, xp] = 0.0
            if xi[xp] < 0.0 or xi[xp] > nx-1 or yi[yp] < 0.0 or yi[yp] > ny-1:
                yi[yp, xp] = np.nan
                continue
            x1 = int(xi[xp]) - a + 1
            if x1 < 0:
                x1 = 0
            if x1 + a2 - 1 > nx:
                x1 = nx - a2 + 1
            x2 = x1 + a2 - 1

            y1 = int(yi[yp]) - a + 1
            if y1 < 0:
                y1 = 0
            if y1 + a2 - 1 > ny:
                y1 = ny - a2 + 1
            y2 = y1 + a2 - 1
            for x in range(x1, x2):
                if np.abs(xi[xp] - x) < a:
                    wx = np.sinc((xi[xp] - x)/a)
                else:
                    wx = 0.0
                for y in range(y1, y2):
                    if np.abs(yi[yp] - y) < a:
                        wy = np.sinc((yi[yp] - y)/a)
                    else:
                        wy = 0.0
                    out[yp, xp] += array[y, x] * np.sinc((xi[xp] - x)) * wx * np.sinc((yi[yp] - y)) * wy
    return out


def interpol_lin1d(x, y, xout):
    """
    Performs a linear interpolation of an irregularily sampled 1D signal

    :author: Andreas Reigber
    :param x: The abscissa values of the signal
    :type x: 1-D ndarray float
    :param y: The ordinate values of the signal
    :type y: 1D ndarray
    :param xout: The abscissa values where the interpolates are desired
    :type xout: 1D ndarray float
    :returns: The interpolated values
    """
    f = sp.interpolate.interp1d(x, y, kind='linear')
    return f(xout)


def interpol_lin2d(x, y, z, xout, yout):
    """
    Performs a linear interpolation of a 2D matrix/image irregularily sampled on both axes

    :author: Andreas Reigber
    :param x: The x values of the first axis of z
    :type x: 1-D ndarray float
    :param y: The y values of the second axis of z
    :type y: 1-D ndarray float
    :param z: The input matrix
    :type z: 2D ndarray
    :param xout: The values on the first axis where the interpolates are desired
    :type xout: 1D ndarray float
    :param yout: The values on the second axis where the interpolates are desired
    :type yout: 1D ndarray float
    :returns: The interpolated matrix /  image
    """
    if np.iscomplexobj(z):
        f_real = sp.interpolate.interp2d(x, y, z.real, kind='linear')
        f_imag = sp.interpolate.interp2d(x, y, z.imag, kind='linear')
        return f_real(xout, yout) + 1j * f_imag(xout, yout)
    else:
        f = sp.interpolate.interp2d(x, y, z, kind='linear')
        return f(xout, yout)


try:
    from .interpolation_extensions import cinterpol_cubic as interpol_cubic
    from .interpolation_extensions import cinterpol_cubic_irr as interpol_cubic_irr
    from .interpolation_extensions import cinterpol_lanczos as interpol_lanczos
    from .interpolation_extensions import cinterpol2D_lanczos as interpol2D_lanczos
except ImportError:
    print("STEtools: Fast cython interpolators not found. (run build process?)")


