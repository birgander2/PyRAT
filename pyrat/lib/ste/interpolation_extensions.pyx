# distutils: extra_compile_args = -fopenmp
# distutils: extra_link_args = -fopenmp

import cython
cimport cython
import numpy as np
cimport numpy as np
from cython.parallel import prange
from libc.math cimport sin, M_PI

ctypedef fused myfloat:
    np.float32_t
    np.float64_t

ctypedef fused fltcpl_t:
    cython.float
    cython.double
    cython.floatcomplex
    cython.doublecomplex

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cinterpol_cubic(fltcpl_t [:] y, myfloat[:] xi, int threads=0):

    cdef int i
    cdef int n = len(y)
    if cython.float is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='complex128')

    cdef fltcpl_t [:] yi = foo
    cdef float nan = np.nan

    cdef fltcpl_t a, b, c, d
    cdef float c00, c10, c20, c30, t, t2, t3
    cdef int khi, klo, kmi, kpl
    cdef int lxi = len(xi)

    for i in prange(lxi, nogil=True, num_threads=threads):
        if xi[i] < 0.0 or xi[i] > n-1:
            yi[i] = nan
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
    return np.asarray(yi)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cinterpol_cubic_irr(myfloat[:] x, fltcpl_t [:] y, myfloat[:] xi, sort=True, int threads=0):
    cdef int i
    cdef int n = len(y)

    if cython.float is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='complex128')

    cdef fltcpl_t [:] yi = foo
    cdef float nan = np.nan

    cdef fltcpl_t a, b, c, d
    cdef float c00, c10, c20, c30, t, t2, t3, h, h2
    cdef int k, khi, klo, kmi, kpl
    cdef int lxi = len(xi)

    if sort is True:
        sidx = np.argsort(x)
        x = np.asarray(x)[sidx]
        y = np.asarray(y)[sidx]

    for i in prange(lxi, nogil=True, num_threads=threads):
        if xi[i] < x[0] or xi[i] > x[n-1]:
            yi[i] = nan
        else:
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
                yi[i] = nan

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
    return np.asarray(yi)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cinterpol_lanczos(fltcpl_t [:] y, myfloat[:] xi, int a=2, int threads=0):
    cdef int i
    cdef int n = len(y)
    if cython.float is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo = np.empty((xi.shape[0], ), dtype='complex128')
    cdef fltcpl_t [:] yi = foo
    cdef int a2, x1, x2, x, klo
    cdef float w, dx
    cdef float nan = np.nan
    cdef float pi = np.pi
    cdef int lxi = len(xi)

    if a < 1:
        a = 1
    a2 = 2*a + 1
    if n < a2:
        print("ERROR: window size too large")
        return
    for i in prange(lxi, nogil=True, num_threads=threads):
        yi[i] = 0.0
        if xi[i] < 0.0 or xi[i] > n-1:
            yi[i] = nan
            continue
        klo = int(xi[i])
        x1 = klo - a + 1
        if x1 < 0:
            x1 = 0
        if x1 + a2 - 1 > n:
            x1 = n - a2 + 1
        x2 = x1 + a2 - 1
        for x in range(x1, x2):
            dx = xi[i] - x
            if dx == 0:
                w = 1.0
            elif dx > a or dx < -a:
                w = 0.0
            else:
                w = sin(pi*dx) * sin(pi*dx/a)*a/pi/pi/dx/dx
            yi[i] = yi[i] + y[x] * w
    return np.asarray(yi)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cinterpol2D_lanczos(fltcpl_t [:, :] array, myfloat[:] yi, myfloat[:] xi, int a=2, int threads=0):
    cdef int ny = array.shape[0]
    cdef int nx = array.shape[1]
    if cython.float is fltcpl_t:
        foo = np.empty((len(yi), len(xi)), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo = np.empty((len(yi), len(xi)), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo = np.empty((len(yi), len(xi)), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo = np.empty((len(yi), len(xi)), dtype='complex128')
    cdef fltcpl_t [:, :] out = foo
    cdef int a2, x1, x2, y1, y2, x, y, klo, xp=0, yp=0
    cdef float wx, wy, dx, dy
    cdef float nan = np.nan
    cdef float pi = np.pi
    cdef int lxi = len(xi)
    cdef int lyi = len(yi)
    if a < 1:
        a = 1
    a2 = 2*a + 1
    if ny < a or nx < a:
        print("ERROR: window size too large")
        return

    for xp in prange(lxi, nogil=True, num_threads=threads):
        if xi[xp] < 0.0 or xi[xp] > nx-1:
            out[:, xp] = nan
            continue
        klo = int(xi[xp])
        x1 = klo - a + 1
        if x1 < 0:
            x1 = 0
        if x1 + a2 - 1 > nx:
            x1 = nx - a2 + 1
        x2 = x1 + a2 - 1
        for yp in range(lyi):
            if yi[yp] < 0.0 or yi[yp] > ny-1:
                out[yp, xp] = nan
                continue
            klo = int(yi[yp])
            y1 = klo - a + 1
            if y1 < 0:
                y1 = 0
            if y1 + a2 - 1 > ny:
                y1 = ny - a2 + 1
            y2 = y1 + a2 - 1
            for x in range(x1, x2):
                dx = xi[xp] - x
                if dx == 0:
                    wx = 1.0
                elif  dx > a or dx < -a:
                    wx = 0.0
                else:
                    wx = sin(pi*dx) * sin(pi*dx/a)*a/pi/pi/dx/dx
                for y in range(y1, y2):
                    dy = yi[xp] - y
                    if dy == 0:
                        wy = 1.0
                    elif  dy > a or dy < -a:
                        wy = 0.0
                    else:
                        wy = sin(pi*dy) * sin(pi*dy/a)*a/pi/pi/dy/dy
                    out[yp, xp] = out[yp, xp] + array[y, x] * wx * wy
    return out


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cinterpol2D_cubic(fltcpl_t [:, :] arr, myfloat[:] yi, myfloat[:] xi, int threads=0):

    cdef int k, l
    cdef int nry = arr.shape[0]
    cdef int nrx = arr.shape[1]

    if cython.float is fltcpl_t:
        foo1 = np.empty((nry, xi.shape[0]), dtype='float32')
        foo2 = np.empty((yi.shape[0], xi.shape[0]), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo1 = np.empty((nry, xi.shape[0]), dtype='complex64')
        foo2 = np.empty((yi.shape[0], xi.shape[0]), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo1 = np.empty((nry, xi.shape[0]), dtype='float64')
        foo2 = np.empty((yi.shape[0], xi.shape[0]), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo1 = np.empty((nry, xi.shape[0]), dtype='complex128')
        foo2 = np.empty((yi.shape[0], xi.shape[0]), dtype='complex128')

    cdef fltcpl_t [:, :] arri1 = foo1
    cdef fltcpl_t [:, :] arri2 = foo2
    cdef float nan = np.nan

    cdef fltcpl_t a, b, c, d
    cdef float c00, c10, c20, c30, t, t2, t3
    cdef int khi, klo, kmi, kpl
    cdef int lxi = len(xi)
    cdef int lyi = len(yi)

    for k in prange(nry, nogil=True, num_threads=threads):
        for l in range(lxi):
            if xi[l] < 0.0 or xi[l] > nrx-1:
                arri1[k, l] = nan
            else:
                klo = int(xi[l])
                khi = klo + 1

                if float(klo) == xi[l]:
                    arri1[k, l] = arr[k, klo]
                else:
                    kmi = klo - 1
                    kpl = khi + 1

                    if kmi < 0:
                        a = 3 * arr[k, klo] - 3 * arr[k, khi] + arr[k, kpl]
                    else:
                        a = arr[k, kmi]

                    if kpl > nrx-1:
                        d = 3 * arr[k, khi] - 3 * arr[k, klo] + arr[k, kmi]
                    else:
                        d = arr[k, kpl]

                    b = arr[k, klo]
                    c = arr[k, khi]
                    t = (xi[l] - klo)
                    t2 = t * t
                    t3 = t2 * t
                    c00 = (- t3 + 2 * t2 - t) / 2.0
                    c10 = (3 * t3 - 5 * t2 + 2) / 2.0
                    c20 = (- 3 * t3 + 4 * t2 + t) / 2.0
                    c30 = (t3 - t2) / 2.0
                    arri1[k, l] = a * c00 + b * c10 + c * c20 + d * c30

    for l in prange(lxi, nogil=True, num_threads=threads):
        for k in range(lyi):
            if yi[k] < 0.0 or yi[k] > nry-1:
                arri2[k, l] = nan
            else:
                klo = int(yi[k])
                khi = klo + 1

                if float(klo) == yi[k]:
                    arri2[k, l] = arri1[klo, l]
                else:
                    kmi = klo - 1
                    kpl = khi + 1

                    if kmi < 0:
                        a = 3 * arri1[klo, l] - 3 * arri1[khi, l] + arri1[kpl, l]
                    else:
                        a = arri1[kmi, l]

                    if kpl > nry-1:
                        d = 3 * arri1[khi, l] - 3 * arri1[klo, l] + arri1[kmi, l]
                    else:
                        d = arri1[kpl, l]

                    b = arri1[klo, l]
                    c = arri1[khi, l]
                    t = (yi[k] - klo)
                    t2 = t * t
                    t3 = t2 * t
                    c00 = (- t3 + 2 * t2 - t) / 2.0
                    c10 = (3 * t3 - 5 * t2 + 2) / 2.0
                    c20 = (- 3 * t3 + 4 * t2 + t) / 2.0
                    c30 = (t3 - t2) / 2.0
                    arri2[k, l] = a * c00 + b * c10 + c * c20 + d * c30

    return np.asarray(arri2)


