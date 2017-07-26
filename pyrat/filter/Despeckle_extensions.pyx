import cython
cimport cython

import scipy as sp

import numpy as np
cimport numpy as np
from libc.math cimport sqrt, abs, exp, lgamma, tgamma
from scipy.ndimage.filters import median_filter



ctypedef fused fltcpl_t:
    cython.float
    cython.double
    cython.floatcomplex
    cython.doublecomplex

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_bilateral(float [:, :] array, win=11, looks=1.0):
    cdef float [:, :] out = np.zeros_like(array)
    cdef float sigma_d = sqrt(win*win/4.6051)
    cdef float sigma_r = np.sqrt(1.0/looks)
    
    cdef int ny = array.shape[0]
    cdef int nx = array.shape[1]
    cdef int ym = win/2
    cdef int xm = win/2
    
    cdef int k, l, x, y
    cdef float norm = 0.0
    cdef float dd, dr, sd, sr
    for k in range(ym, ny-ym):
        for l in range(xm, nx-xm):
            norm = 0.0
            for y in range(-ym, ym+1):
                for x in range(-xm, xm+1):
                    dd = sqrt(x**2 + y**2)                     # spatial distance
                    dr = abs(array[k, l] - array[k+y, l+x])    # statistical distance
                    sd = exp(-(dd / sigma_d)**2 / 2)
                    sr = exp(-(dr / sigma_r)**2 / 2)
                    out[k, l] += array[k+y, l+x] * sd * sr
                    norm += sd * sr
            out[k, l] /= norm
    return np.asarray(out)

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_leesigmaold(float [:, :] array, float looks=1.0, win=(7,7)):
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
def cy_leesigma(float [:, :] span, fltcpl_t [:, :, :, :] array, float looks=1.0, win=(7,7)):
    cdef fltcpl_t [:, :, :, :] out = np.empty_like(array)
    cdef float diff = 2.0 * np.sqrt(1.0/looks)
    cdef int ny = array.shape[2]
    cdef int nx = array.shape[3]
    cdef int nv = array.shape[0]
    cdef int nz = array.shape[1]
    cdef int ym = win[0]/2
    cdef int xm = win[1]/2
    cdef int limit = (xm + ym)  # //4
    cdef int k, l, x, y, v, z
    if cython.float is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='complex128')
    cdef fltcpl_t [:, :] res = foo

    cdef int n = 0
    for k in range(ym, ny-ym):
        for l in range(xm, nx-xm):
            for v in range(nv):
                for z in range(nz):
                    res[v, z] = 0.0
            n = 0
            for y in range(-ym, ym+1):
                for x in range(-xm, xm+1):
                    if span[k+y, l+x]>span[k, l]*(1.0-diff) and (span[k+y, l+x]<span[k, l]*(1.0+diff)):
                        for v in range(nv):
                            for z in range(nz):
                                res[v, z] = res[v, z] + array[v, z, k+y, l+x]
                        n += 1
            if n >= limit:
                for v in range(nv):
                    for z in range(nz):
                        out[v, z, k, l] = res[v, z] / n
            else:
               for v in range(nv):
                    for z in range(nz):
                        out[v, z, k, l] = (array[v, z, k-1, l] + array[v, z, k+1, l] + array[v, z, k, l-1] + array[v, z, k, l+1]) / 4.0
    return np.asarray(out)


# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_leeimproved_old(float [:, :] array, bounds=(0.5, 3.0), float thres=5.0, looks=1.0, win=(9, 9), float newsig=0.5):
    cdef float sig2 = 1.0 / looks
    cdef float sfak = 1.0 + sig2
    cdef float nsig2 = newsig
    cdef float nsfak = 1.0 + nsig2
    cdef int ny = array.shape[0]
    cdef int nx = array.shape[1]
    cdef int ym = win[0]/2
    cdef int xm = win[1]/2
    cdef int k, l, x, y
    cdef float m2arr, marr, vary, varx, kfac, i1, i2
    cdef float res = 0.0
    cdef int n = 0

    cdef float [:, :] out = np.zeros_like(array)
    for k in range(ym, ny-ym):
        for l in range(xm, nx-xm):
            m2arr = 0.0
            marr = 0.0
            n = 0
            for y in range(-1, 2):                          # check 3x3 neighbourhood
                for x in range(-1, 2):
                    m2arr += array[k+y, l+x]**2
                    marr += array[k+y, l+x]
                    if array[k+y, l+x] > thres:
                        n += 1

            if n >= 6:                                      # keep all point targets
                for y in range(-1, 2):
                    for x in range(-1, 2):
                        if array[k+y, l+x] > thres:
                            out[k+y, l+x] = array[k+y, l+x]

            if out[k, l] == 0.0:                             # no point target
                m2arr /= 9.0
                marr /= 9.0
                vary = (m2arr - marr**2)
                if vary < 1e-10: vary = 1e-10
                varx = ((vary - marr ** 2 * sig2) / sfak)
                if varx < 0: varx = 0
                kfac = varx / vary
                xtilde = (array[k, l] - marr) * kfac + marr

                i1 = xtilde*bounds[0]
                i2 = xtilde*bounds[1]
                m2arr = 0.0
                marr = 0.0
                n = 0

                for y in range(-ym, ym+1):
                    for x in range(-xm, xm+1):
                        if array[k+y, l+x]>i1 and array[k+y, l+x]<i2:
                            m2arr += array[k+y, l+x]**2
                            marr += array[k+y, l+x]
                            n += 1
                if n == 0:
                    out[k, l] = 0.0
                else:
                    m2arr /= n
                    marr /= n
                    vary = (m2arr - marr**2)
                    if vary < 1e-10: vary = 1e-10
                    varx = ((vary - marr ** 2 * nsig2) / nsfak)
                    if varx < 0.0: varx = 0.0
                    kfac = varx / vary
                    out[k, l] = (array[k, l] - marr) * kfac + marr
    return np.asarray(out)

# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_leeimproved(float [:, :] span, fltcpl_t [:, :, :, :] array, bounds=(0.5, 3.0), float thres=5.0, looks=1.0, win=(9, 9), float newsig=0.5):
    cdef fltcpl_t [:, :, :, :] out = np.zeros_like(array)
    cdef float sig2 = 1.0 / looks
    cdef float sfak = 1.0 + sig2
    cdef float nsig2 = newsig
    cdef float nsfak = 1.0 + nsig2
    cdef float xtilde
    cdef int nv = array.shape[0]
    cdef int nz = array.shape[1]
    cdef int ny = array.shape[2]
    cdef int nx = array.shape[3]
    cdef int ym = win[0]/2
    cdef int xm = win[1]/2
    cdef int norm = win[0] * win[1]
    cdef int k, l, x, y, v, z
    cdef float m2arr, marr, vary, varx, kfac, i1, i2
    if cython.float is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='complex128')
    cdef fltcpl_t [:, :] res = foo

    cdef int n = 0
    for k in range(ym, ny-ym):
        for l in range(xm, nx-xm):
            m2arr = 0.0
            marr = 0.0
            n = 0
            for y in range(-1, 2):                          # check 3x3 neighbourhood
                for x in range(-1, 2):
                    m2arr += span[k+y, l+x]**2
                    marr += span[k+y, l+x]
                    if span[k+y, l+x] > thres:
                        n += 1
            if n >= 6:                                      # keep all point targets
                for y in range(-1, 2):
                    for x in range(-1, 2):
                        for v in range(nv):
                            for z in range(nz):
                                if span[k+y, l+x] > thres:
                                    out[v, z, k+y, l+x] = array[v, z, k+y, l+x]

            if out[0, 0, k, l] == 0.0:                      # no point target, also not prior
                m2arr /= 9.0
                marr /= 9.0
                vary = (m2arr - marr**2)
                if vary < 1e-10: vary = 1e-10
                varx = ((vary - marr ** 2 * sig2) / sfak)
                if varx < 0: varx = 0
                kfac = varx / vary

                xtilde = (span[k, l] - marr) * kfac + marr

                i1 = xtilde*bounds[0]
                i2 = xtilde*bounds[1]
                m2arr = 0.0
                marr = 0.0
                n = 0
                for v in range(nv):
                    for z in range(nz):
                        res[v, z] = 0.0

                for y in range(-ym, ym+1):
                    for x in range(-xm, xm+1):
                        if span[k+y, l+x]>i1 and span[k+y, l+x]<i2:
                            m2arr += span[k+y, l+x]**2
                            marr += span[k+y, l+x]
                            n += 1
                            for v in range(nv):
                                for z in range(nz):
                                    res[v, z] = res[v, z] + array[v, z, k+y, l+x]
                if n == 0:
                    for v in range(nv):
                        for z in range(nz):
                            out[v, z, k, l] = 0.0
                else:
                    m2arr /= n
                    marr /= n
                    vary = (m2arr - marr**2)
                    if vary < 1e-10: vary = 1e-10
                    varx = ((vary - marr ** 2 * nsig2) / nsfak)
                    if varx < 0.0: varx = 0.0
                    kfac = varx / vary
                    for v in range(nv):
                        for z in range(nz):
                            out[v, z, k, l] = (array[v, z, k, l] - res[v, z] / n) * kfac + res[v, z] / n
    return np.asarray(out)


# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_srad_old(float [:, :] array, float looks=1.0, float step=0.05, int iter=0, float scale=1.0):
    cdef float p1, p2, p3, p4, aim, aip, ajm, ajp
    cdef int i, j
    cdef float d2i, qp1, qp2, q2, d
    cdef int ny = array.shape[0]
    cdef int nx = array.shape[1]

    cdef q0 = 1.0/sqrt(looks)*exp(-step*iter/6.0)

    cdef float [:, :] ci = np.zeros_like(array)
    cdef float [:, :] out = np.zeros_like(array)

    for i in range(0, ny):
        for j in range(0, nx):
            aip = array[i+1, j] if i != ny-1 else array[i, j]
            aim = array[i-1, j] if i != 0 else array[i, j]
            ajp = array[i, j+1] if j != nx-1 else array[i, j]
            ajm = array[i, j-1] if j != 0 else array[i, j]
            p1 = aip - array[i, j]
            p2 = ajp - array[i, j]
            p3 = array[i, j] - aim
            p4 = array[i, j] - ajm

            d2I = (aip + aim + ajp + ajm - 4.0*array[i, j]) / scale / scale
            qp1 = sqrt(p1*p1 + p2*p2 + p3*p3 + p4*p4) / array[i, j] / scale
            qp2 = d2I / array[i, j]

            q = sqrt((qp1*qp1/2.0 - qp2*qp2/16.0) / (1.0 + qp2/4.0)**2)

            ci[i, j] = 1/(1+(q*q - q0*q0)/(q0*q0*(1.0+q0*q0)))
            # ci[i, j] = exp(-(q*q - q0*q0)/(q0*q0*(1.0+q0*q0)))
            # ci[i, j] = 1.0

    for i in range(0, ny):
        for j in range(0, nx):
            aip = array[i+1, j] if i != ny-1 else array[i, j]
            aim = array[i-1, j] if i != 0 else array[i, j]
            ajp = array[i, j+1] if j != nx-1 else array[i, j]
            ajm = array[i, j-1] if j != 0 else array[i, j]
            cip = ci[i+1, j] if i != ny-1 else ci[i, j]
            cim = ci[i-1, j] if i != 0 else ci[i, j]
            cjp = ci[i, j+1] if j != nx-1 else ci[i, j]
            cjm = ci[i, j-1] if j != 0 else ci[i, j]

            d = (cip*(aip-array[i, j]) + ci[i, j]*(aim-array[i, j]) + cjp*(ajp-array[i, j]) + ci[i, j]*(ajm-array[i, j]))
            out[i, j] = array[i, j] + step/4.0*d

    return np.asarray(out)


# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_srad(float [:, :] span, fltcpl_t [:, :, :, :] array, float looks=1.0, float step=0.05, int iter=0, float scale=1.0):
    cdef float p1, p2, p3, p4, sip, sim, sjp, sjm
    cdef fltcpl_t aim, aip, ajm, ajp
    cdef int i, j, v, z
    cdef float d2i, qp1, qp2, q2
    cdef int nv = array.shape[0]
    cdef int nz = array.shape[1]
    cdef int ny = array.shape[2]
    cdef int nx = array.shape[3]

    cdef q0 = 1.0/sqrt(looks)*exp(-step*iter/6.0)
    cdef q02 = q0 * q0

    cdef float [:, :] ci = np.zeros_like(span, dtype='f4')
    cdef fltcpl_t [:, :, :, :] out = np.zeros_like(array)

    for i in range(1, ny-1):
        for j in range(1, nx-1):
            p1 = span[i+1, j] - span[i, j]
            p2 = span[i, j+1] - span[i, j]
            p3 = span[i, j] - span[i-1, j]
            p4 = span[i, j] - span[i, j-1]
            d2I = (span[i+1, j] + span[i-1, j] + span[i, j+1] + span[i, j-1] - 4.0*span[i, j]) / scale / scale
            qp1 = sqrt(p1**2 + p2**2 + p3**2 + p4**2) / span[i, j] / scale
            qp2 = d2I / span[i, j]
            q = ((qp1*qp1/2.0 - qp2*qp2/16.0) / (1.0 + qp2/4.0)**2)
            ci[i, j] = 1/(1+(q - q02)/(q02*(1.0+q02)))
            # ci[i, j] = exp(-(q - q02)/(q02*(1.0+q02)))
            # ci[i, j] = 1.0

    # boundary conditions 1

    for i in [0, ny-1]:
        for j in [0, nx-1]:
            sip = span[i+1, j] if i != ny-1 else span[i, j]
            sim = span[i-1, j] if i != 0 else span[i, j]
            sjp = span[i, j+1] if j != nx-1 else span[i, j]
            sjm = span[i, j-1] if j != 0 else span[i, j]
            p1 = sip - span[i, j]
            p2 = sjp - span[i, j]
            p3 = span[i, j] - sim
            p4 = span[i, j] - sjm
            d2I = (sip + sim + sjp + sjm - 4.0*span[i, j]) / scale / scale
            qp1 = sqrt(p1**2 + p2**2 + p3**2 + p4**2) / span[i, j] / scale
            qp2 = d2I / span[i, j]
            q = ((qp1*qp1/2.0 - qp2*qp2/16.0) / (1.0 + qp2/4.0)**2)
            ci[i, j] = 1/(1+(q - q02)/(q02*(1.0+q02)))
            # ci[i, j] = exp(-(q - q02)/(q02*(1.0+q02)))
            # ci[i, j] = 1.0

    # update

    for i in range(1, ny-1):
        for j in range(1, nx-1):
            for v in range(nv):
                for z in range(nz):
                    out[v, z, i, j] = array[v, z, i, j] + step/4.0*(
                          ci[i+1, j]*(array[v, z, i+1, j]-array[v, z, i, j])
                        + ci[i, j]*(array[v, z, i-1, j]-array[v, z, i, j])
                        + ci[i, j+1]*(array[v, z, i, j+1]-array[v, z, i, j])
                        + ci[i, j]*(array[v, z, i, j-1]-array[v, z, i, j]))

    # boundary conditions 2

    for i in [0, ny-1]:
        for j in [0, nx-1]:
            cip = ci[i+1, j] if i != ny-1 else ci[i, j]
            cjp = ci[i, j+1] if j != nx-1 else ci[i, j]
            for v in range(nv):
                for z in range(nz):
                    aip = array[v, z, i+1, j] if i != ny-1 else array[v, z, i, j]
                    aim = array[v, z, i-1, j] if i != 0 else array[v, z, i, j]
                    ajp = array[v, z, i, j+1] if j != nx-1 else array[v, z, i, j]
                    ajm = array[v, z, i, j-1] if j != 0 else array[v, z, i, j]
                    out[v, z, i, j] = array[v, z, i, j] + step/4.0*(cip*(aip-array[v, z, i, j])
                        + ci[i, j]*(aim-array[v, z, i, j]) + cjp*(ajp-array[v, z, i, j])
                        + ci[i, j]*(ajm-array[v, z, i, j]))

    return np.asarray(out)

# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_emdes(float [:, :] span, fltcpl_t [:, :, :, :] array, float looks=1.0, win=(7,7)):
    cdef fltcpl_t [:, :, :, :] out = np.empty_like(array)
    cdef float sig2 = 1.0 / looks
    cdef float sfak = 1.0 + sig2
    cdef int ny = array.shape[2]
    cdef int nx = array.shape[3]
    cdef int nv = array.shape[0]
    cdef int nz = array.shape[1]
    cdef int ym = win[0]/2
    cdef int xm = win[1]/2
    cdef int limit = (xm + ym)  # //4
    cdef int k, l, x, y, v, z
    if cython.float is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='complex128')
    cdef fltcpl_t [:, :] res = foo
    cdef float p, pi, fl, p1, p2

    for k in range(ym, ny-ym):
        for l in range(xm, nx-xm):

            # --------------- estimate centre pixel intensity ---------------------------

            m2arr = 0.0
            marr = 0.0
            for y in range(-1, 2):                          # check 3x3 neighbourhood
                for x in range(-1, 2):
                    m2arr += span[k+y, l+x]**2
                    marr += span[k+y, l+x]
            m2arr /= 9.0
            marr /= 9.0
            vary = (m2arr - marr**2)
            if vary < 1e-10: vary = 1e-10
            varx = ((vary - marr ** 2 * sig2) / sfak)
            if varx < 0: varx = 0
            kfac = varx / vary
            xtilde = (span[k, l] - marr) * kfac + marr

            xtilde = span[k, l]
            # --------------- estimate centre pixel intensity ---------------------------
            for v in range(nv):
                for z in range(nz):
                    res[v, z] = 0.0
            pi = 0.0

            for y in range(-ym, ym+1):
                for x in range(-xm, xm+1):

#                    fl = looks - (span[k+y, l+x] - xtilde)/(span[k+y, l+x] + xtilde)
#                    if fl <= 1.0:
#                        fl = 1.0
#                    p1 = (span[k+y, l+x]**fl * xtilde**fl)
#                    p2 = ((span[k+y, l+x] + xtilde)/2)**(2*fl)
#                    if p2 != 0.0:
#                        p = (p1 / p2)
#                    else:
#                        p = 0.0

                    if xtilde != 0.0:
                        p = (looks**looks * xtilde**(looks-1.0))/(tgamma(looks)*span[k+y, l+x]**looks)*exp(-looks*xtilde/span[k+y, l+x])
                        # p = (looks**looks * span[k+y, l+x]**(looks-1.0))/(tgamma(looks)*xtilde**looks)*exp(-looks*span[k+y, l+x]/xtilde)
                    else:
                        p = 0.0

                    pi  += p
                    for v in range(nv):
                        for z in range(nz):
                            res[v, z] = res[v, z] + p * array[v, z, k+y, l+x]

            for v in range(nv):
                for z in range(nz):
                    if pi != 0.0:
                        out[v, z, k, l] = res[v, z] / pi

    return np.asarray(out)

# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def cy_idanq(span, fltcpl_t [:, :, :, :] array, float looks=1.0, int nmax=50, bint llmmse=True):
    cdef fltcpl_t [:, :, :, :] out = np.empty_like(array)
    cdef float sig2 = 1.0 / looks
    cdef float sfak = 1.0 + sig2
    cdef float klee, varx, vary, imean
    cdef float limit = 1.0 / sqrt(looks) * 2.0 / 3.0
    cdef int ny = array.shape[2]
    cdef int nx = array.shape[3]
    cdef int nv = array.shape[0]
    cdef int nz = array.shape[1]
    cdef float [:] seeds = median_filter(span, (3, 3)).flatten()

    span[0, :] = 1e10
    span[-1, :] = 1e10
    span[:, 0] = 1e10
    span[:, -1] = 1e10
    cdef float [:] amp = span.flatten()

    if cython.float is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='float32')
    elif cython.floatcomplex is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='complex64')
    elif cython.double is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='float64')
    elif cython.doublecomplex is fltcpl_t:
        foo = np.zeros((nv, nz), dtype='complex128')
    cdef fltcpl_t [:, :] res = foo
    cdef int nb3[8]
    cdef unsigned long region[10000]
    cdef unsigned long backgr[10000]
    cdef int i, j, x, y, k, l
    cdef int nold, nnew, nbak
    cdef unsigned long ix, idx, pos
    cdef float seed, limit3 = limit*3.0
    cdef bint check
    nb3[:] = [-nx-1, -nx, -nx+1, -1, 1, nx+1, nx, nx-1]
    for x in range(1, nx-1):
        for y in range(1, ny-1):
            pos = y * nx + x
            seed = seeds[pos]
            region[0] = pos
            nold = 1
            nbak = 0

            while True:
                nnew = nold
                for i in range(nold):
                    idx = region[i]
                    for j in range(8):
                        ix = idx + nb3[j]
                        check = False
                        for k in range(nnew):
                            if ix == region[k]:
                                check = True
                        if check is False:
                            if (abs(amp[ix] - seed) / seed) < limit:
                                region[nnew] = ix
                                nnew += 1
                            else:
                                check = False
                                for k in range(nbak):
                                    if ix == backgr[k]:
                                        check = True
                                if check is False:
                                    backgr[nbak] = ix
                                    nbak += 1
                if nnew >= nmax:
                    break
                if nnew == nold:
                    break
                nold = nnew

            # update seed
            seed = 0.0
            for ix in region[0:nnew]:
                seed += amp[ix]
            seed /= nnew

            # check remaining background pixels
            for ix in backgr[0:nbak]:
                if (abs(amp[ix] - seed) / seed) < limit3:
                    region[nnew] = ix
                    nnew += 1

            for v in range(nv):
                for z in range(nz):
                    res[v, z] = 0.0

            for v in range(nv):
                for z in range(nz):
                    for ix in region[0:nnew]:
                        res[v, z] = res[v, z] + array[v, z, ix // nx, ix % nx]

            if llmmse == True:
                # calculate LLMMSE (Lee) factor
                imean = 0.0
                vary = 0
                for ix in region[0:nnew]:
                    imean += amp[ix]
                imean /= nnew
                for ix in region[0:nnew]:
                    vary += (amp[ix] - imean)**2
                vary /= (nnew - 1)
                varx = ((vary - sig2 * imean ** 2) / sfak)
                if varx < 0.0:
                    varx = 0.0
                klee = varx / vary

                # LLMMSE filtering of region
                for v in range(nv):
                    for z in range(nz):
                        out[v, z, y, x] = (array[v, z, y, x] - res[v, z] / nnew) * klee + res[v, z] / nnew
            else:
                for v in range(nv):
                    for z in range(nz):
                        out[v, z, y, x] = res[v, z] / nnew

    return np.asarray(out)

