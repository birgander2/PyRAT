import PyRat, glob, os, logging
import numpy as np
import ipdb

class ESAR_track(PyRat.ImportWorker):
    """
    Import of DLR E-SAR REF/REALTRACK products
    
    :author: Andreas Reigber
    :param filename: The filename of the DCSLC file to import.
    :type array: string
    :param setlen: Interpolate to the given azimuth length.
    :type array: int
  
    """
    def __init__(self, *args, **kwargs):
        super(ESAR_track, self).__init__(*args, **kwargs)    
        self.name  = "ESAR REAL/REF TRACK IMPORT"
        self.block = 'T'
        
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
                nraz = self.setlen
                track= np.empty((nraz,7),dtype='float64')
                time = np.interp(np.arange(nraz,dtype='float64')/nraz*nry,np.arange(nry,dtype='float64'),time) 
                for l in range(3):
                    track[:,l+0] = np.interp(np.arange(nraz,dtype='float64')/nraz*nry,np.arange(nry,dtype='float64'),trk2[l,:]) 
                    track[:,l+3] = np.interp(np.arange(nraz,dtype='float64')/nraz*nry,np.arange(nry,dtype='float64'),trk1[l,:]) 
                track[:,6] = time
            else:
                track= np.empty((nry,7),dtype='float64')
                for l in range(3):
                    track[:,l+0] = trk2[l,:]
                    track[:,l+3] = trk1[l,:]
                track[:,6] = time
            return track, meta
        else:
            logging.error("Track file not found (setting track to None): "+self.filename)
            return None, None
            
