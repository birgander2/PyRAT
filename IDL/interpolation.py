import numpy as np
from scipy.ndimage import interpolation
from scipy import interpolate as sci_int
from scipy import ndimage
import scipy.spatial.qhull as qhull

import IDL


def shift2(arr, *args):
    """
    Like normal (IDL) shift, but allows subpixel offsets. Uses shift when offsets are integer, interpolates when float.                                                                                                                           

    :author: Andreas Reigber
    :param arr: The input array (numpy ndarray)
    :param offsets: list of offsets, same lengths as number of dimensions of the input array
    :type arr1: int or float
    :returns: shifted input array
    """
    size = arr.shape    
    if np.equal(np.floor(args),np.array(args)).all():
        oarr = IDL.shift(arr, *args)
    else:
        if np.iscomplexobj(arr):
            oarr = interpolation.shift(arr.real,args) + 1j * interpolation.shift(arr.imag,args)
        else:
            oarr = interpolation.shift(arr,args)
    return oarr 




def interpol(y, x, xnew, spline=0):
    # real and imaginary parts separately...
    if (np.iscomplexobj(y)):
        return interpol(y.real,x,xnew,spline=spline) \
            + 1j*interpol(y.imag,x,xnew,spline=spline)

    if (spline):
        spl_coef = sci_int.splrep(x,y,s=0)
        return np.asarray(sci_int.splev(xnew,spl_coef,der=0),dtype=y.dtype)
    else:
        return np.asarray(np.interp(xnew,x,y),dtype=y.dtype)




def interpolate(arr,*args,**kwargs):
    """
    Interpolates a n-dimensional array. The first argument is the array to be
    interpolated. The next arguments are, respectively, the positions in the
    i-th dimension onto which the arr should be interpolated.

    If the argument grid is set to True, then the positions are vectors
    defining a n-dimensional grid.

    If the argument cubic is set to True a cubic interpolation is performed

    The function replicates the behaviour of interpolate in IDL.
    :param arr: the input n-dimensional array
    """

    if (len(args) != arr.ndim):
        raise Exception('Expected the same number of coordinate matrices as array dimensions!')

    if (np.iscomplexobj(arr)):
        return interpolate(arr.real,*args,**kwargs) \
            + 1j*interpolate(arr.imag,*args,**kwargs)

    if ('grid' in kwargs) and (kwargs['grid']):
        # re-gridding onto a regular grid
        arrShape = np.asarray(arr.shape)
        ndim = arr.ndim
        for n in range(ndim):
            view = np.rollaxis(arr,n,ndim)
            lastShape = np.asarray(view.shape)
            view.reshape((view.size/view.shape[ndim-1],view.shape[ndim-1]))
            x = np.arange(view.shape[1])
            xnew = np.asarray(args[n]).flat

            outShape = (view.shape[0], len(xnew))
            if (outShape[1] == view.shape[1] and n > 0):
                # in-place interpolation
                out = view
            else:
                # need new array (in the first iteration always)
                out = np.empty(outShape,dtype=view.dtype)

            if ('cubic' in kwargs) and (kwargs['cubic'] != 0):
                for r in xrange(view.shape[0]):
                    spl_coef = sci_int.splrep(x,view[r,:],s=0)
                    out[r, :] = sci_int.splev(xnew,spl_coef,der=0)
            else:
                for r in xrange(view.shape[0]):
                    out[r, :] = np.interp(xnew,x,view[r,:])

            lastShape[ndim-1] = len(xnew)
            arr = np.rollaxis(out.reshape(lastShape),ndim-1,n)

        return arr
    else:
        order = 2
        if ('cubic' in kwargs) and (kwargs['cubic'] != 0):
            order = 3

        return ndimage.map_coordinates(arr,args,order=order,output=arr.dtype)



def triangulate(xyz, uvw):
    """
    Found on the internet. Not quite the same but works...
    Prepares point cloud(s) for linear interpolation.

    xyz: irregular points [Npoints x Ndims] (source)
    uvw: irregular or regular points [Npoints x Ndims] (destination)
    """
    
    d = xyz.shape[1]
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))



def trigrid(values, (vtx, wts), missing=np.nan):
    """
    Found on the internet. Not quite the same but works...
    Performs linear interpolation on 'values' to re-grid from xyz -> uvw (see triangulate above)

    values: the scalar value at xyz (see triangulate above)
    (vtx,wts): output of triangulate
    missing: value to substitute for samples outside the convex hull of xyz
    """

    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = missing
    return ret


