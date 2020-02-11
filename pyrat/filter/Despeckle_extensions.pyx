# cython: language_level=3
import cython
cimport cython

import scipy as sp

import numpy as np
cimport numpy as np
from cpython cimport array
from libc.math cimport sqrt, abs, exp, lgamma, tgamma
from scipy.ndimage.filters import median_filter

ctypedef fused fltcpl_t:
    cython.float
    cython.double
    cython.floatcomplex
    cython.doublecomplex

ctypedef cython.doublecomplex DTYPE_t

cdef extern from "math.h":
#    double sqrt(double)
    double sin(double)
    double cos(double)
    double atan2(double, double)
    double fabs(double)
    double log(double)
#    double exp(double)
cdef extern from "complex.h":
    double creal(double complex)
    double cimag(double complex)
    double complex conj(double complex)
    double complex clog(double complex)
cdef extern from "eigenH.h":
    double SQR(double)
    double SQR_ABS(double complex)
    cdef double DBL_EPSILON "_DBL_EPSILON"
    cdef double sqrt3 "_sqrt3"

# Main Beltrami routine. Calculates the neighbors distances using a region growing algorithm.
def cy_MCB(np.ndarray[DTYPE_t, ndim=3] cov_array, np.ndarray[np.float32_t, ndim=2] distance_array,  np.ndarray[np.int_t,
    ndim=2] all_neigh, np.ndarray[np.int_t, ndim=1] neighbours_local, float sigma, float beta, float betastr,
     int window_size, long rdim, looks, bint llmmse=True):
    cdef int ddim = np.shape(cov_array)[1]
    cdef double complex [:,:,:] avg = np.empty_like(cov_array)
    cdef double fac = sigma/(betastr * beta)
    cdef int maxp = np.shape(cov_array)[0]
    cdef int i, j, k, z, iter1, iter2, dx
    cdef int totalW = window_size * window_size
    cdef int sizeNeigh = np.shape(all_neigh)[0]
    cdef long [:] nx
    cdef int flag
    cdef long center_pixel2
    cdef long localID, local_idx
    cdef double [:] dBeltrami = np.zeros([totalW]) + np.inf
    cdef long [:] idxBeltrami = np.zeros([totalW], dtype=np.int64) - 30
    cdef int middle = (totalW) // 2
    cdef double original_distance
    cdef long coord
    cdef double [:] exp1
    cdef double complex [:,:] up = np.zeros((ddim, ddim), dtype = np.complex_)
    cdef double down
    cdef double sigma2 = sigma**2
    cdef double [:] tmp = np.zeros(totalW)
    cdef double complex [:,:,:] included = np.zeros((totalW, ddim, ddim), dtype = np.complex_)
    cdef float currentD, tempD
    cdef float gamma = 1
    cdef double gammad = np.sqrt(2)
    cdef double klee, varx, vary, imean, vary0
    cdef float sig2 = 1.0 / looks
    cdef float sfak = 1.0 + sig2

    # For all pixels
    for k in range(maxp):
        for i in range(totalW):
            # Initilizating the distances with a big number
            dBeltrami[i] = 100000000
            # Initializing the indexes to avoid edge issues
            idxBeltrami[i] = - 30
        middle = (window_size * window_size) // 2
        dBeltrami[middle], idxBeltrami[middle] = 0, k
        for i in range(totalW):
            tmp[i] = 0

        # For all pixels inside the local window (NxN)
        for iter1 in range(sizeNeigh):
            nx = all_neigh[iter1]
            flag = 0
            center_pixel2 = k + nx[1]
            localID = middle + nx[0]
            original_distance = dBeltrami[localID]
            # For all pixels inside the immediate window (closest 8 pixels). The algorithm obtains a previous (infinite)
            # distance and then increases it depending on the distance to the central pixel and the
            for iter2 in range(8):
                dx = neighbours_local[iter2]
                local_idx = localID + dx
                currentD = dBeltrami[local_idx]
                coord = center_pixel2
                if flag == 0:
                    if 0 <= coord < maxp:
                        tempD = gammad + fac*distance_array[0, coord] + original_distance
                        if tempD < currentD:
                            dBeltrami[local_idx] = tempD
                            if 0 <= coord - rdim + 1 < maxp:
                                idxBeltrami[local_idx] = coord - rdim + 1

                elif flag == 1:
                    if 0 <= coord < maxp:
                        tempD = gamma + fac*distance_array[1, coord] + original_distance
                        if tempD < currentD:
                            dBeltrami[local_idx] = tempD
                            if 0 <= coord + 1 < maxp:
                                idxBeltrami[local_idx] = coord + 1

                elif flag == 2:
                    if 0 <= coord < maxp:
                        tempD = gammad + fac*distance_array[2, coord] + original_distance
                        if tempD < currentD:
                            dBeltrami[local_idx] = tempD
                            if 0 <= coord + 1 + rdim < maxp:
                                idxBeltrami[local_idx] = coord + 1 + rdim

                elif flag == 3:
                    if 0 <= coord < maxp:
                        tempD = gamma + fac*distance_array[3, coord] + original_distance
                        if tempD < currentD:
                            dBeltrami[local_idx] = tempD
                            if 0 <= coord + rdim < maxp:
                                idxBeltrami[local_idx] = coord + rdim

                elif flag == 4:
                    coord = center_pixel2 - 1 + rdim
                    if 0 <= coord < maxp:
                        tempD = gammad + fac*distance_array[0, coord] + original_distance
                        if tempD < currentD:
                            dBeltrami[local_idx] = tempD
                            idxBeltrami[local_idx] = coord
                elif flag == 5:
                    coord = center_pixel2 - 1
                    if 0 <= coord < maxp:
                        tempD = gamma + fac*distance_array[1, coord] + original_distance
                        if tempD < currentD:
                            dBeltrami[local_idx] = tempD
                            idxBeltrami[local_idx] = coord
                elif flag == 6:
                    coord = center_pixel2 - rdim - 1
                    if 0 <= coord < maxp:
                        tempD = gammad + fac*distance_array[2, coord] + original_distance
                        if tempD < currentD:
                            dBeltrami[local_idx] = tempD
                            idxBeltrami[local_idx] = coord
                elif flag == 7:
                    coord = center_pixel2 - rdim
                    if 0 <= coord < maxp:
                        tempD = gamma + fac*distance_array[3, coord] + original_distance
                        if tempD < currentD:
                            dBeltrami[local_idx] = tempD
                            idxBeltrami[local_idx] = coord
                flag += 1
            for i in range(totalW):
                if idxBeltrami[i] == -30:
                    idxBeltrami[i] = center_pixel2

        included = cov_array[idxBeltrami, :, :]
        for i in range(totalW):
            tmp[i] = exp((-dBeltrami[i]**2) / sigma2)
        down = 0
        for i in range(ddim):
            for j in range(ddim):
                 up[i, j] = 0

        for z in range(totalW):
            for i in range(ddim):
                for j in range(ddim):
                    up[i, j] += included[z,i,j] * tmp[z]

        for i in range(totalW):
            down += tmp[i]

        if llmmse == True:

            # calculate LLMMSE (Lee) factor
            imean = 0.0
            vary = 0
            vary0 = 0
            for i in range(totalW):
                for j in range(ddim):
                    imean += creal(included[i,j,j])
            imean /= totalW

            for i in range(totalW):
                for j in range(ddim):
                    vary0 += creal(included[i,j,j])
                vary -= vary0 - (imean)**2
            vary /= (totalW - 1)
            varx = ((vary - sig2 * imean ** 2) / sfak)
            if varx < 0.0:
                varx = 0.0
            klee = varx / vary

            # LLMMSE filtering of region
            for i in range(ddim):
                for j in range(ddim):
                    avg[k, i, j] = (cov_array[k,i,j] - up[i, j] / down) * klee + up[i, j] / down
            #for v in range(nv):
            #    for z in range(nz):
            #        out[v, z, y, x] = (array[v, z, y, x] - res[v, z] / nnew) * klee + res[v, z] / nnew
        else:
            #        out[v, z, y, x] = res[v, z] / nnew
            for i in range(ddim):
                for j in range(ddim):
                    avg[k, i, j] = up[i, j] / down
    return np.squeeze(avg)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
# Cython routine for the fast Affine-invariant distance.
def cy_MCBdist(np.ndarray[DTYPE_t, ndim=2] A, np.ndarray[DTYPE_t, ndim=2] B):
    cdef double complex [3][3] Sig1 = A
    cdef double complex [3][3] Sig2 = B
    cdef double complex [3][3] sqrt_invA
    cdef double complex [3][3] v
    cdef double complex [3][3] tmp
    cdef double complex [3][3] tmp2
    cdef double complex [3][3] m1
    cdef double complex [3][3] m2
    cdef double complex [3][3] m3
    cdef double [3] w
    cdef double out

    cdef int i, j
    #Initializing...
    for i in range(3):
        for j in range(3):
            sqrt_invA[i][j] = 0

    w = eigVal(Sig1)
    v = eigVec(Sig1, w)
    # Obtaining the inverse sqrt of A
    for i in range(3):
        for j in range(3):
            if (j == i): sqrt_invA[i][j] = sqrt(1/w[i])

    tmp = matMulH(sqrt_invA, v)
    m1 = matMul(v,tmp)
    tmp2 = matMul(Sig2, m1)
    m2 = matMul(m1, tmp2)
    m3 = matLog(m2)
    out = frob(m3)
    return out

def cy_eigen(np.ndarray[DTYPE_t, ndim=2] B):
    cdef double complex [3][3] A = B
    cdef double complex [3][3] Q
    cdef double [3] w
    w = eigVal(A)
    Q = eigVec(A, w)
    return np.asarray(w), np.asarray(Q)

cdef matLog(double complex [3][3] a):
    cdef double complex [3][3] diagPrime
    cdef double complex [3][3] tmp2
    cdef double complex [3][3] tmp3
    cdef double complex [3][3] tmp0
    cdef double complex [3][3] logm
    cdef double complex [3][3] out
    cdef double complex [3][3] aprime
    cdef double complex [3][3] v
    cdef double [3] w
    cdef int i, j
    cdef int x, y
    #Initializing...
    for x in range(3):
        for y in range(3):
            diagPrime[x][y] = 0 + 0j

    w = eigVal(a)
    v = eigVec(a, w)
    tmp0 = matMul(a, v)
    aprime = matHMul(v, tmp0)

    for i in range(3):
        for j in range(3):
            if (j == i):
                diagPrime[i][j] = clog(aprime[i][j])
    tmp2 = matMulH(diagPrime, v)
    logm = matMul(v, tmp2)
    return logm


cdef frob(double complex [3][3] a):
    cdef double result = 0
    cdef double complex value = 0
    cdef int i,j
    for i in range(3):
        for j in range(3):
            value = a[i][j]
            result += SQR_ABS(value)
    return sqrt(result)

cdef matMul(double complex [3][3] a, double complex [3][3] b):
    #Simple matrix multiplication
    cdef double complex [3][3] c
    cdef int i, j, k
    cdef int x, y
    #Initializing...
    for x in range(3):
        for y in range(3):
            c[x][y] = 0 + 0j

    for i in range(3):
        for j in range(3):
            for k in range(3):
                c[i][j] = c[i][j] + (a[i][k]*b[k][j])
    return c

cdef matMulH(double complex [3][3] a, double complex [3][3] b):
    #Matrix multiplication with conj(b).T as the second input
    cdef double complex [3][3] c
    cdef int i, j, k

    #Initializing...
    for i in range(3):
        for j in range(3):
            c[i][j] = 0 + 0j
    for i in range(3):
        for j in range(3):
            for k in range(3):
                c[i][j] += (a[i][k]*conj(b[j][k]))
    return c

cdef matHMul(double complex [3][3] a, double complex [3][3] b):
    #Matrix multiplication with conj(b).T as the second input
    cdef double complex [3][3] c
    cdef int i, j, k

    #Initializing...
    for i in range(3):
        for j in range(3):
            c[i][j] = 0 + 0j

    for i in range(3):
        for j in range(3):
            for k in range(3):
                c[i][j] += (conj(a[k][i])*(b[k][j]))
    return c

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef eigVal(double complex [3][3] B):
    cdef double complex [3][3] A = B
    cdef double [3] w
    cdef double m, c1, c0, p, sqrt_p, q, c, s, phi
    cdef double complex de = A[0][1] * A[1][2]
    cdef double dd = SQR_ABS(A[0][1])
    cdef double ee = SQR_ABS(A[1][2])
    cdef double ff = SQR_ABS(A[0][2])

    m  = creal(A[0][0]) + creal(A[1][1]) + creal(A[2][2]);
    c1 = (creal(A[0][0])*creal(A[1][1])  + creal(A[0][0])*creal(A[2][2]) + creal(A[1][1])*creal(A[2][2]))- (dd + ee + ff)
    c0 = creal(A[2][2])*dd + creal(A[0][0])*ee + creal(A[1][1])*ff - creal(A[0][0])*creal(A[1][1])*creal(A[2][2])- 2.0 * (creal(A[0][2])*creal(de) + cimag(A[0][2])*cimag(de))

    p = SQR(m) - 3.0*c1
    q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0
    sqrt_p = sqrt(fabs(p))

    phi = 27.0 * ( 0.25*SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0))
    phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q)

    c = sqrt_p*cos(phi);
    s = (1.0/sqrt3)*sqrt_p*sin(phi)

    w[1]  = (1.0/3.0)*(m - c)
    w[2]  = w[1] + s
    w[0]  = w[1] + c
    w[1] -= s
    return w

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef eigVec(double complex [3][3] B, double [3] w):
    cdef double complex [3][3] A = B
    cdef double complex [3][3] Q
    cdef double norm, n0, n1, n0tmp, n1tmp, thresh, error, wmax
    cdef int i, j
    cdef double complex f

    wmax = fabs(w[0])
    if (fabs(w[1]) > wmax):
        wmax = fabs(w[1])
    if (fabs(w[2]) > wmax):
        wmax = fabs(w[2])

    thresh = SQR(8.0 * DBL_EPSILON * wmax)

    # Prepare calculation of eigenvectors
    n0tmp   = SQR_ABS(A[0][1]) + SQR_ABS(A[0][2])
    n1tmp   = SQR_ABS(A[0][1]) + SQR_ABS(A[1][2])
    Q[0][1] = A[0][1]*A[1][2] - A[0][2]*creal(A[1][1])
    Q[1][1] = A[0][2]*conj(A[0][1]) - A[1][2]*creal(A[0][0])
    Q[2][1] = SQR_ABS(A[0][1])

    #Calculate first eigenvector by the formula
    # v[0] = conj( (A - w[0]).e1 x (A - w[0]).e2 )

    A[0][0] = A[0][0] - w[0]
    A[1][1] = A[1][1] - w[0]
    Q[0][0] = Q[0][1] + A[0][2]*w[0]
    Q[1][0] = Q[1][1] + A[1][2]*w[0]
    Q[2][0] = creal(A[0][0])*creal(A[1][1]) - Q[2][1]
    norm    = SQR_ABS(Q[0][0]) + SQR_ABS(Q[1][0]) + SQR(creal(Q[2][0]))
    n0      = n0tmp + SQR(creal(A[0][0]))
    n1      = n1tmp + SQR(creal(A[1][1]))
    error   = n0 * n1

    if (n0 <= thresh):  # If the first column is zero, then (1,0,0) is an eigenvector
        Q[0][0] = 1.0
        Q[1][0] = 0.0
        Q[2][0] = 0.0
    elif (n1 <= thresh):# If the second column is zero, then (0,1,0) is an eigenvector
        Q[0][0] = 0.0
        Q[1][0] = 1.0
        Q[2][0] = 0.0
    elif (norm < SQR(64.0 * DBL_EPSILON) * error):
                                # If angle between A[0] and A[1] is too small, don't use
        t = SQR_ABS(A[0][1])   # cross product, but calculate v ~ (1, -A0/A1, 0)
        f = -A[0][0] / A[0][1]
        if (SQR_ABS(A[1][1]) > t):
          t = SQR_ABS(A[1][1])
          f = -conj(A[0][1]) / A[1][1]
        if (SQR_ABS(A[1][2]) > t):
          f = -conj(A[0][2] / A[1][2])
        norm = 1.0/sqrt(1 + SQR_ABS(f))
        Q[0][0] = norm
        Q[1][0] = f * norm
        Q[2][0] = 0.0
    else:
        norm = sqrt(1.0 / norm)
        for j in range(3):
            Q[j][0] = Q[j][0] * norm

#    Prepare calculation of second eigenvector
    t = w[0] - w[1]
    if (fabs(t) > 8.0 * DBL_EPSILON * wmax):
        # For non-degenerate eigenvalue, calculate second eigenvector by the formula
        #   v[1] = conj( (A - w[1]).e1 x (A - w[1]).e2 )
        A[0][0]  = A[0][0] + t
        A[1][1]  = A[1][1] + t
        Q[0][1]  = Q[0][1] + A[0][2]*w[1]
        Q[1][1]  = Q[1][1] + A[1][2]*w[1]
        Q[2][1]  = creal(A[0][0])*creal(A[1][1]) - creal(Q[2][1])
        norm     = SQR_ABS(Q[0][1]) + SQR_ABS(Q[1][1]) + SQR(creal(Q[2][1]))
        n0       = n0tmp + SQR(creal(A[0][0]))
        n1       = n1tmp + SQR(creal(A[1][1]))
        error    = n0 * n1

        if (n0 <= thresh):       # If the first column is zero, then (1,0,0) is an eigenvector
          Q[0][1] = 1.0
          Q[1][1] = 0.0
          Q[2][1] = 0.0
        elif (n1 <= thresh):  # If the second column is zero, then (0,1,0) is an eigenvector
          Q[0][1] = 0.0
          Q[1][1] = 1.0
          Q[2][1] = 0.0
        elif (norm < SQR(64.0 * DBL_EPSILON) * error):
                                   # If angle between A[0] and A[1] is too small, don't use
          t = SQR_ABS(A[0][1])   # cross product, but calculate v ~ (1, -A0/A1, 0)
          f = -A[0][0] / A[0][1]
          if (SQR_ABS(A[1][1]) > t):
            t = SQR_ABS(A[1][1])
            f = -conj(A[0][1]) / A[1][1]
          if (SQR_ABS(A[1][2]) > t):
            f = -conj(A[0][2] / A[1][2])
            norm = 1.0/sqrt(1 + SQR_ABS(f))

          Q[0][1] = norm
          Q[1][1] = f * norm
          Q[2][1] = 0.0
        else:
          norm = sqrt(1.0 / norm)
          for j in range(3):
              Q[j][1] = Q[j][1] * norm
    else:
      # For degenerate eigenvalue, calculate second eigenvector according to
      #   v[1] = conj( v[0] x (A - w[1]).e[i] )

      # This would really get to complicated if we could not assume all of A to
      # contain meaningful values.
      A[1][0]  = conj(A[0][1])
      A[2][0]  = conj(A[0][2])
      A[2][1]  = conj(A[1][2])
      A[0][0]  = A[0][0] + w[0]
      A[1][1]  = A[1][1] + w[0]
      for i in range(3):
         A[i][i] = A[i][i] - w[1]
         n0        = SQR_ABS(A[0][i]) + SQR_ABS(A[1][i]) + SQR_ABS(A[2][i])
         if (n0 > thresh):
            Q[0][1]  = conj(Q[1][0]*A[2][i] - Q[2][0]*A[1][i])
            Q[1][1]  = conj(Q[2][0]*A[0][i] - Q[0][0]*A[2][i])
            Q[2][1]  = conj(Q[0][0]*A[1][i] - Q[1][0]*A[0][i])
            norm     = SQR_ABS(Q[0][1]) + SQR_ABS(Q[1][1]) + SQR_ABS(Q[2][1])
            if (norm > SQR(256.0 * DBL_EPSILON) * n0):  # Accept cross product only if the angle between
                                                        # the two vectors was not too small
              norm = sqrt(1.0 / norm)

              for j in range(3):
                Q[j][1] = Q[j][1] * norm

    # Calculate third eigenvector according to
    #   v[2] = conj(v[0] x v[1])
    Q[0][2] = conj(Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1])
    Q[1][2] = conj(Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1])
    Q[2][2] = conj(Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1])

    return Q


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

