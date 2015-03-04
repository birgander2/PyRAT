import pyrat, glob, os, logging
import numpy as np


class ESAR_track(pyrat.ImportWorker):
    """
    Import of DLR E-SAR REF/REALTRACK products
    
    :author: Andreas Reigber
    :param filename: The filename of the DCSLC file to import.
    :type array: string
    :param setlen: Interpolate to the given azimuth length.
    :type array: int
    :param crop: Sets a crop area to read. Default: (0,0,0,0) = entire file
    :type crop: tuple
  
    """
    def __init__(self, *args, **kwargs):
        super(ESAR_track, self).__init__(*args, **kwargs)    
        self.name  = "ESAR REAL/REF TRACK IMPORT"
        self.block = 'T'
        if 'crop' not in self.__dict__:  self.crop = (0,0,0,0)
        
    def reader(self, *args, **kwargs):
        meta  = {}
        meta['sensor'] = 'DLR E-SAR'
        array = False
        
        if os.path.isfile(self.filename):
            lun = open(self.filename,'rb')
            nry  = np.fromfile(file=lun, dtype="int32",count=1).byteswap()[0]    
            time = np.fromfile(file=lun, dtype="float64", count=nry).byteswap()
            trk1 = np.fromfile(file=lun, dtype="float64", count=nry*3).byteswap().reshape(3, nry)
            trk2 = np.fromfile(file=lun, dtype="float64", count=nry*3).byteswap().reshape(3, nry)
            trk3 = np.fromfile(file=lun, dtype="float64", count=nry*3).byteswap().reshape(3, nry)
            lun.close()
            
            if hasattr(self, 'setlen'):
                naz = self.setlen
                track= np.empty((naz,7),dtype='float64')
                time = np.interp(np.arange(naz,dtype='float64')/naz*nry,np.arange(nry,dtype='float64'),time) 
                for l in range(3):
                    track[:,l+0] = np.interp(np.arange(naz,dtype='float64')/naz*nry,np.arange(nry,dtype='float64'),trk2[l,:]) 
                    track[:,l+3] = np.interp(np.arange(naz,dtype='float64')/naz*nry,np.arange(nry,dtype='float64'),trk1[l,:]) 
                track[:,6] = time
            else:
                naz = nry
                track= np.empty((naz,7),dtype='float64')
                for l in range(3):
                    track[:,l+0] = trk2[l,:]
                    track[:,l+3] = trk1[l,:]
                track[:,6] = time
            
            block  = list(self.crop)
            if block[1] == 0: block[1] = naz
            return track[block[0]:block[1],:], meta
        else:
            logging.error("Track file not found (setting track to None): "+self.filename)
            return None, None
            
def esar_track(*args, **kwargs):
    return ESAR_track(*args, **kwargs).run(**kwargs)
