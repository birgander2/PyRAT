import PyRat
import numpy as np
import STEtools as STE
import pdb, logging
from scipy.interpolate import interp2d
from scipy import ndimage

class Coregister(PyRat.FilterWorker):
    """
    Lexicographic to Pauli conversion...
    """

    def __init__(self, *args, **kwargs):
        super(Coregister, self).__init__(*args, **kwargs)    
        self.name = "PATCH COREGISTRATION"
        
    def filter(self, array, *args, **kwargs):
        arr1 = array[0]
        arr2 = array[1]
        ny, nx = arr1.shape
        dy, dx = 2**(int(np.log(min(ny, 4096))/np.log(2))),2**(int(np.log(min(nx, 4096))/np.log(2)))
        offset = STE.coreg(arr1[(ny-dy)/2:(ny+dy)/2, (nx-dx)/2:(nx+dx)/2],arr2[(ny-dy)/2:(ny+dy)/2, (nx-dx)/2:(nx+dx)/2])
        logging.info('Global offset : '+str(int(offset[0]))+' / '+str(int(offset[1])))
        
        paty = np.arange(12)*ny/12
        patx = np.arange(12)*nx/12
        dy   = 2**(int(np.log(ny/12)/np.log(2)))
        dx   = 2**(int(np.log(nx/12)/np.log(2)))
        offy = np.zeros((10,10))
        offx = np.zeros((10,10))
        
        for y, yp in enumerate(paty[1:11]):
            for x, xp in enumerate(patx[1:11]):
                amp1 = np.abs(arr1[yp:yp+dy, xp:xp+dx])
                amp2 = np.abs(arr2[yp-int(offset[0]):yp+dy-int(offset[0]), xp-int(offset[1]):xp+dx-int(offset[1])])
                #pdb.set)trace
                foo  = STE.coreg(amp1,amp2, sub=True)
                offy[y,x] = foo[0] + offset[0]
                offx[y,x] = foo[1] + offset[1]
        
        #pdb.set_trace()       
        xx, yy = np.meshgrid(patx[1:11], paty[1:11])
        cx = polyfit2d(yy.flatten(), xx.flatten(), offx.flatten(), order=2)
        cy = polyfit2d(yy.flatten(), xx.flatten(), offy.flatten(), order=2)
        #cx = polyfit2d(xx.flatten(), yy.flatten(), offx.flatten(), order=3)
        #cy = polyfit2d(xx.flatten(), yy.flatten(), offy.flatten(), order=3)
        
        ny, nx = arr2.shape
        xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
        px = polyval2d(yy, xx, cx)
        py = polyval2d(yy, xx, cy)
        #px = polyval2d(xx, yy, cx)
        #py = polyval2d(xx, yy, cy)
        
        arr2 = ndimage.map_coordinates(arr2, np.rollaxis(np.dstack([yy-py,xx-px]),2))
        #pdb.set_trace()
        #arr2 = np.roll(np.roll(arr2, int(offset[0]), axis=0), int(offset[1]), axis=1)
            
        return arr2
    
import itertools
def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    #x = x * 1.0
    #y = y * 1.0
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x, dtype='float64')
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z        