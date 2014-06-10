import PyRat
import scipy as sp
from scipy.signal import convolve2d
from scipy.ndimage import filters
from STEtools import *

import pdb, logging, time
import numpy as np

class RefinedLee(PyRat.FilterWorker):
    """
    Refined Lee speckle filter
    
    further information:
    J.S.Lee et al.: "Speckle Reduction in Multipolarization Multifrequency SAR Imagery", 
    Trans. on Geoscience and Remote Sensing, Vol. 29, No. 4, pp. 535-544, 1991
    """

    def __init__(self, *args, **kwargs):
        super(RefinedLee, self).__init__(*args, **kwargs)    
        self.name = "REFINED LEE FILTER"

        if 'win'    not in self.__dict__: self.win=7
        if 'looks'  not in self.__dict__: self.looks=1.0
        if 'threshold'  not in self.__dict__: self.threshold=0.5
        if 'method' not in self.__dict__: self.method='original'
        self.blockprocess = True
        #self.blocksize    = 512
        #self.nthreads     = 1
        self.blockoverlap = self.win/2+1
        
    def filter(self, array, *args, **kwargs):

# ---------------------------------------------
# INIT & SPAN
# ---------------------------------------------

        sig2  = 1.0 / self.looks
        sfak  = 1.0 + sig2
        nrx   = array.shape[-1]

        lshape = array.shape[0:-2]
        if len(lshape) == 2:
            #span = np.abs(np.trace(array,axis1=0,axis2=1))
            span = np.abs(array[0,0,...] + array[1,1,...] + array[2,2,...])
        else:
            logging.error("Data not in matrix form")
    
# ---------------------------------------------
# TURNING BOX 
# ---------------------------------------------

        cbox = np.zeros((9,self.win,self.win),dtype='float32')
        chbox = np.zeros((self.win,self.win),dtype='float32')
        chbox[0:self.win//2+1,:] = 1
        cvbox = np.zeros((self.win,self.win),dtype='float32')
        for k in range(self.win):
            cvbox[k,0:k+1] = 1
        
        cbox[0,...] = np.rot90(chbox,3)
        cbox[1,...] = np.rot90(cvbox,1)
        cbox[2,...] = np.rot90(chbox,2)
        cbox[3,...] = np.rot90(cvbox,0)
        cbox[4,...] = np.rot90(chbox,1)
        cbox[5,...] = np.rot90(cvbox,3)
        cbox[6,...] = np.rot90(chbox,0)
        cbox[7,...] = np.rot90(cvbox,2)
        for k in range(self.win//2+1):
            for l in range(self.win//2-k,self.win//2+k+1):
                cbox[8,k:self.win-k,l] = 1

        for k in range(9): 
            cbox[k,...] /= np.sum(cbox[k,...])
        
        ampf1 = np.empty((9,)+span.shape)
        ampf2 = np.empty((9,)+span.shape)
        for k in range(9):
            ampf1[k,...] = filters.correlate(span**2,cbox[k,...])
            ampf2[k,...] = filters.correlate(span,cbox[k,...])**2
    
# ---------------------------------------------
# GRADIENT ESTIMATION
# ---------------------------------------------
        np.seterr(divide='ignore',invalid='ignore')
        if self.method == 'original':
            xs = [+2,+2, 0,-2,-2,-2, 0,+2]
            ys = [0 ,+2,+2,+2, 0,-2,-2,-2]
            samp = filters.uniform_filter(span,self.win//2)
            grad = np.empty((8,)+span.shape)
            for k in range(8):
                grad[k,...] = np.abs(np.roll(np.roll(samp, ys[k], axis=0), xs[k], axis=1)/samp - 1.0)
            magni = np.amax(grad,axis=0)
            direc = np.argmax(grad,axis=0)
            direc[magni < self.threshold] = 8
        elif self.method == 'cov':
            grad = np.empty((8,)+span.shape)
            for k in range(8):
                grad[k,...] = np.abs((ampf1[k,...] - ampf2[k,...]) / ampf2[k,...])
                direc = np.argmin(grad,axis=0)
        else:
            logging.error("Unknown method!")  
        
        np.seterr(divide='warn',invalid='warn')
# ---------------------------------------------
# FILTERING
# ---------------------------------------------
        out  = np.empty_like(array)
        dbox = np.zeros((1,1)+(self.win,self.win))
        for l in range(9):
            grad = ampf1[l,...]
            mamp = ampf2[l,...]
            dbox[0,0,...] = cbox[l,...]
            
            vary = (grad - mamp).clip(1e-10)
            varx = ((vary - mamp*sig2)/sfak).clip(0)
            kfac = varx / vary
            mamp = filters.correlate(array.real,dbox) + 1j*filters.convolve(array.imag,dbox)

            idx = np.argwhere(direc == l)            
            out[:,:,idx[:,0],idx[:,1]] = (mamp + (array-mamp)*kfac)[:,:,idx[:,0],idx[:,1]]
        
        return out
        
        
        
        
        
        