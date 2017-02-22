# distutils: extra_compile_args = -fopenmp
# distutils: extra_link_args = -fopenmp

import numpy as np
cimport numpy as np
ctypedef np.float32_t dtype_t
import cython
cimport cython
from cython.parallel import prange
from libc.math cimport sqrt, complex

@cython.boundscheck(False)
@cython.wraparound(False)
def mtrebin2d(float [:, :] arr, float [:, :] out, int threads=0):
    cdef int by = out.shape[0]
    cdef int bx = out.shape[1]
    cdef int ny = arr.shape[0]
    cdef int nx = arr.shape[1]
    cdef int fy = arr.shape[0] / by
    cdef int fx = arr.shape[1] / bx
    cdef float norm = fy * fx
    cdef int y=0, x=0, dy=0, dx=0

    for y in prange(by, nogil=True, num_threads=threads):
        for x in range(bx):
            for dy in range(fy):
                for dx in range(fx):
                    out[y,x] += arr[y*fy+dy, x*fx+dx]
            out[y, x] /= norm

def parallel_rebin(arr, size, **kwargs):
    dshape = tuple(arr.shape[0:2])
    lshape = tuple(arr.shape[2:])
    nchannels = int(np.prod(lshape))
    if nchannels == 1:
        arr = arr[:, :, np.newaxis]
    else:
        arr = arr.reshape(dshape+(nchannels, ))
        size = size[0:2]
    out = np.zeros(tuple(size) + (nchannels,), dtype=arr.dtype)
    for k in range(nchannels):
        mtrebin2d(arr[..., k], out[..., k], **kwargs)
    return out.reshape(tuple(size)+lshape)

@cython.boundscheck(False)
@cython.wraparound(False)
def mtabs2d(float complex [:, :] arr, int threads=8):
    cdef int ny = arr.shape[0]
    cdef int nx = arr.shape[1]
    cdef float [:, :] out = np.empty((ny, nx), dtype='f4')
    cdef int y=0, x=0

    for y in prange(ny, nogil=True, num_threads=threads):
        for x in range(nx):
            out[y,x] = sqrt(arr[y, x].real**2 + arr[y, x].imag**2)
    return np.asarray(out)

def parallel_abs(arr, **kwargs):
    shp = arr.shape
    dshape = tuple(shp[0:2])
    lshape = tuple(shp[2:])
    nchannels = int(np.prod(lshape))
    if nchannels == 1:
        arr = arr[:, :, np.newaxis]
    else:
        arr = arr.reshape(dshape+(nchannels, ))
    out = np.empty_like(arr, dtype=np.float32)
    for k in range(nchannels):
        out[..., k] = mtabs2d(arr[..., k], **kwargs)
    return out.reshape(shp)

