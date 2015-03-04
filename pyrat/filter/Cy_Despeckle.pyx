import cython
cimport cython

import scipy as sp

import numpy as np
cimport numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_leesigma(float [:, :] array, float looks=1.0, win=(7,7)):
    cdef float [:, :] out = np.empty_like(array)
    cdef float diff = 2.0 * np.sqrt(1.0/looks)
    cdef int ny = array.shape[0]
    cdef int nx = array.shape[1]
    cdef int ym = win[0]/2
    cdef int xm = win[1]/2
    cdef int limit = (xm + ym)  # //4
    cdef int k, l, x, y
    cdef float res = 0.0
    cdef int n = 0
    for k in range(ym, ny-ym):
        for l in range(xm, nx-xm):
            res = 0.0
            n = 0
            for y in range(-ym, ym+1):
                for x in range(-xm, xm+1):
                    if array[k+y, l+x]>array[k, l]*(1.0-diff) and (array[k+y, l+x]<array[k, l]*(1.0+diff)):
                        res += array[k+y, l+x]
                        n += 1
            if n >= limit:
                out[k, l] = res / n
            else:
                out[k, l] = (array[k-1, l] + array[k+1, l] + array[k, l-1] + array[k, l+1])/4.0
    return np.asarray(out)


@cython.boundscheck(False)
@cython.wraparound(False)
def cy_leesigmanew(float [:, :] array, float looks=1.0, win=(7,7)):

# STEP 2

    sig2 = 1.0 / looks
    sfak = 1.0 + sig2
    span = np.asarray(array)
    m2arr = sp.ndimage.filters.uniform_filter(span**2, size=(5,5))
    marr = sp.ndimage.filters.uniform_filter(span, size=(5,5))
    vary = (m2arr - marr ** 2).clip(1e-10)
    varx = ((vary - marr ** 2 * sig2) / sfak).clip(0)
    kfac = varx / vary
    xestnp = marr + (span - marr) * kfac

# STEP 3
    cdef float diff = 2.0 * np.sqrt(1.0/looks)

    cdef float [:, :] out = np.empty_like(array)
    cdef float [:, :] xest = xestnp
    cdef int ny = array.shape[0]
    cdef int nx = array.shape[1]
    cdef int ym = win[0]/2
    cdef int xm = win[1]/2
    cdef int k, l, x, y

    cdef float m2pix = 0.0
    cdef float mpix = 0.0
    cdef int npix = 0
    cdef float k2, varzz, varxx
    cdef i1 = 0.221
    cdef i2 = 2.744
    cdef nv2 = 0.5699**2

    for k in range(ym, ny-ym):
        for l in range(xm, nx-xm):
            npix = 0
            mpix = 0.0
            m2pix = 0.0
            for y in range(-ym, ym+1):
                for x in range(-xm, xm+1):
                    if (array[k+y, l+x]>=xest[k, l]*i1) and (array[k+y, l+x]<=xest[k, l]*i2):
                    # if array[k+y, l+x]>array[k, l]*(1.0-diff) and (array[k+y, l+x]<array[k, l]*(1.0+diff)):
                        mpix += array[k+y, l+x]
                        m2pix += array[k+y, l+x]*array[k+y, l+x]
                        npix += 1
            if npix>0:
                mpix /= npix
                m2pix /= npix
                varzz = m2pix - mpix*mpix
                varxx = (varzz - mpix*mpix*nv2) / (1.0+nv2)
                if varzz>0.0:
                    k2 = varxx / varzz
                else:
                    k2 = 0.0
            else:
                k2 = 1.0
            out[k, l] = mpix  + (array[k, l] - mpix) * k2

    return np.asarray(out)


