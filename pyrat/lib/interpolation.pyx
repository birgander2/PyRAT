# distutils: extra_compile_args = -fopenmp
# distutils: extra_link_args = -fopenmp

import cython
cimport cython
import numpy as np
cimport numpy as np
from cython.parallel import prange

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
def interpol_cubic(fltcpl_t [:] y, myfloat[:] xi, int threads=0):

    cdef int i
    cdef int n = len(y)
    cdef float[:] x = np.arange(n, dtype=np.float32)
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
        if xi[i] < x[0] or xi[i] > x[-1]:
            yi[i] = nan

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
def interpol_cubic_irr(myfloat[:] x, fltcpl_t [:] y, myfloat[:] xi, sort=True, int threads=0):
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
        if xi[i] < x[0] or xi[i] > x[-1]:
            yi[i] = nan

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

