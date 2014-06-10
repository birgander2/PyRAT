import PyRat, glob, os, logging
import numpy as np

class ESAR_dcslc(PyRat.ImportWorker):
    """
    Import of DLR E-SAR DCSLC products
    
    This module imports E-SAR DCSLC files. If a wildcard for the channel number is used,
    the complete polarimetric scattering vector is imported. All files (i.e. SLC data,
    efile and track) must be in the same directory.
    
    :param filename: The filename of the DCSLC file to import.
    :type array: string
    :author: Andreas Reigber
   
    """
    def __init__(self, *args, **kwargs):
        super(ESAR_dcslc, self).__init__(*args, **kwargs)    
        self.name = "ESAR DCSLC IMPORT"
        
    def reader(self, *args, **kwargs):
        meta  = {}
        meta['sensor'] = 'DLR E-SAR'
        array = False
        
        files = glob.glob(self.filename)
        
        if os.path.isfile(files[0]):
            lun = open(files[0],'rb')
            header = np.fromfile(file=lun, dtype="int32",count=3).byteswap()
            lun.close()
            npix = header[1] * header[2]
            nraz = header[1]
            array = np.empty((len(files),header[1],header[2]),dtype='complex64')
            meta['CH_pol'] = [' ']*len(files)
        else:
            logging.error("File not found: "+file)
            return 
        
        for k, file in enumerate(files):
            logging.info("Found "+file)
            
            # READ DATA
            if os.path.isfile(file):
                lun = open(file,'rb')
                header = np.fromfile(file=lun, dtype="int32",count=3).byteswap()
                array[k,...] = np.fromfile(file=lun, dtype="complex64", count=npix).byteswap().reshape(header[1],header[2])
                lun.close()
            else:
                logging.error("File not found: "+file)
            
            # READ EFILE - METADATA
            
            fl = os.path.split(file)
            efile = fl[0] + '/e' + fl[1][1:].replace('_dcslc.dat', '.txt')
            if os.path.isfile(efile):
                lun = open(efile)
                etext = lun.readlines()
                lun.close()
                meta['CH_pol'][k] = next(line for line in etext if 'init.polarization' in line)[-3:-1]
                if k==0:
                    meta['prf'] = float(next(line for line in etext if 'init.prf' in line)[-15:])
                    meta['c0']  = float(next(line for line in etext if 'init.speed_of_light' in line)[-15:])
                    meta['rd']  = float(next(line for line in etext if 'init.range_delay' in line)[-15:])*1e-6
                    meta['rs']  = float(next(line for line in etext if 'init.range_sampling_rate' in line)[-15:])*1e6
                    meta['lambda']  = float(next(line for line in etext if 'init.wavelength' in line)[-15:])
                    meta['band'] = next(line for line in etext if 'init.freq_band' in line)[-15:].strip()
                    meta['antdir'] = -1
            else:
                logging.error("Metadata file not found: "+efile)
        
        return array, meta
