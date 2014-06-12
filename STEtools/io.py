from scipy import misc
import __builtin__
import numpy
import ctypes
import os, string


#-------------------------------------------------------------------------      
def write_pixmap(filename, image_array, palette=False):
    """
    Saves a numpy 2-D ndarray as jpg/png/etc.
    
    Parameters
    ----------
    filename : str
        Output filename.
    image : ndarray, MxN or MxNx3 or MxNx4
        Array containing image values. If the shape is ``MxN``, the array
        represents a grey-level image. Shape ``MxNx3`` stores the red, green
        and blue bands along the last dimension. An alpha layer may be
        included, specified as the last colour band of an ``MxNx4`` array.
    """
    if palette == True and hasattr(__builtin__,'myPalette'):
        p = __builtin__.myPalette
        misc.imsave(filename, p[image_array])
    else:
        misc.imsave(filename, image_array)
        
write_png = write_pixmap
write_jpg = write_pixmap

        
#-------------------------------------------------------------------------      
#-------------------------------------------------------------------------      
# RAT FILE HANDLING
#-------------------------------------------------------------------------      
#-------------------------------------------------------------------------      
def rrat(filename, **kwargs):
    """
    Reads an entire RAT file, returns it as a numpy ndarray variable
    """
    File = RatFile(filename)
    return File.read(**kwargs)
    
#-------------------------------------------------------------------------      
def srat(filename,array,**kwargs):
    """
    Writes a numpy ndarray into a RAT file
    """
    File = RatFile(filename)
    File.write(array,**kwargs)

#-------------------------------------------------------------------------      
class RatHeaderRat(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("magiclong",  ctypes.c_int ),
                ("version",    ctypes.c_float),
                ("ndim",       ctypes.c_int),
                ("nchannel",   ctypes.c_int),
                ("dim",        ctypes.c_int * 8),
                ("var",        ctypes.c_int),
                ("sub",        ctypes.c_int * 2),
                ("rattype",    ctypes.c_int ),
                ("reserved",   ctypes.c_int * 9)]
    def __init__(self):
        self.magiclong = 844382546
        self.version   = 2.0
        self.sub       = (ctypes.c_int * 2)(1, 1)
     
class RatHeaderInfo(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("info",   ctypes.c_char * 100)]

class RatHeaderGeo(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("projection",   ctypes.c_short),
                ("ps_east",      ctypes.c_double),
                ("ps_north",     ctypes.c_double),
                ("min_east",     ctypes.c_double),
                ("min_north",    ctypes.c_double),
                ("zone",         ctypes.c_short),
                ("hemisphere",   ctypes.c_short),
                ("long0scl",     ctypes.c_double),
                ("max_axis_ell", ctypes.c_double),
                ("min_axis_ell", ctypes.c_double),
                ("dshift_tx",    ctypes.c_double),
                ("dshift_ty",    ctypes.c_double),
                ("dshift_tz",    ctypes.c_double),
                ("dshift_rx",    ctypes.c_double),
                ("dshift_ry",    ctypes.c_double),
                ("dshift_rz",    ctypes.c_double),
                ("dshift_scl",   ctypes.c_double),
                ("dshift_info",  ctypes.c_char * 64),
                ("reserved",     ctypes.c_byte * 18)]

class RatHeaderEmpty(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("reserved",   ctypes.c_int * 25)]

class RatHeader(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("Rat",       RatHeaderRat),
                ("Info",      RatHeaderInfo),
                ("Geo",       RatHeaderGeo),
                ("Stat",      RatHeaderEmpty),
                ("Reserved1", RatHeaderEmpty),
                ("Reserved2", RatHeaderEmpty),
                ("Reserved3", RatHeaderEmpty),
                ("Reserved4", RatHeaderEmpty),
                ("Reserved5", RatHeaderEmpty)]
                
    def __init__(self):
        self.Rat = RatHeaderRat()

#-------------------------------------------------------------------------      
#-------------------------------------------------------------------------      
class RatFile():
    """
    Class for manipulating RAT formatted files. Work in progress, many features still missing...
    """
   
    def __init__(self,filename):
        self.filename = filename
        self.Header = RatHeader()
        if os.path.exists(self.filename):
            self.version, self.xdrflag = self.version()
            self.read_header()
            self.exists = True
        else:
            self.exists = False         
        
#--------------------------------------------------------------------------           
    def read(self, **kwargs):
        if self.version == 2.0:
            lun = open(self.filename,'rb')
            if 'block' in kwargs:
                block = kwargs['block']
                blockshape = self.dim[:self.ndim]
                blockshape[1] = block[3]
                arr = numpy.ndarray(shape=blockshape[::-1],dtype=self.dtype)
                lun.seek(1000 + arr.itemsize * block[1] * blockshape[0] * self.nchannel)
                lun.readinto(arr.data)
                arr = arr[:,block[0]:block[0]+block[2]]
            else:
                lun = open(self.filename,'rb')
                lun.seek(1000)
                arr = numpy.ndarray(shape=self.dim[0:self.ndim][::-1],dtype=self.dtype)
                lun.readinto(arr.data)
            lun.close()
            return arr
        elif self.version == 1.0:
            lun = open(self.filename,'rb')
            lun.seek(104+4*self.ndim+4*self.xdrflag)
            arr = numpy.ndarray(shape=self.dim[::-1],dtype=self.dtype)
            lun.readinto(arr.data)
            lun.close()
            if self.xdrflag == 1: arr = arr.byteswap()
            return arr    
        else:
            print "ERROR: RAT version not supported"
            raise IOError

#--------------------------------------------------------------------------   
    def write(self,arr,type=0,**kwargs):
        if ('header' in kwargs):
            self.Header = kwargs['header']

        self.dtype = arr.dtype
        if   arr.dtype == 'uint8':     self.var = 1
        elif arr.dtype == 'int16':     self.var = 2
        elif arr.dtype == 'int32':     self.var = 3
        elif arr.dtype == 'float32':   self.var = 4
        elif arr.dtype == 'float64':   self.var = 5
        elif arr.dtype == 'complex64': self.var = 6
        elif arr.dtype == 'complex128':self.var = 9
        elif arr.dtype == 'uint16':    self.var = 12
        elif arr.dtype == 'uint32':    self.var = 13
        elif arr.dtype == 'int64':     self.var = 14
        elif arr.dtype == 'uint64':    self.var = 15
        self.dim   = arr.shape[::-1]
        self.ndim  = arr.ndim
        self.Header.Rat.ndim = self.ndim
        self.Header.Rat.dim  = self.dim
        self.Header.Rat.var  = self.var
        self.Header.Rat.rattype = type
        self.info = string.rstrip(self.Header.Info.info)
        
        lun = open(self.filename,'wb')
        lun.write(self.Header)
        #lun.write(arr.tostring())
        arr.tofile(lun)
        lun.close()
        #self.help()

#--------------------------------------------------------------------------   
    def read_header(self):      
        if self.version == 2.0:
            lun = open(self.filename,'rb')
            lun.readinto(self.Header)
            lun.close()
            self.ndim     = int(self.Header.Rat.ndim)
            self.nchannel = int(self.Header.Rat.nchannel)
            self.dim      = [int(foo) for foo in self.Header.Rat.dim]
            self.var      = int(self.Header.Rat.var)
            self.sub      = [int(foo) for foo in self.Header.Rat.sub]
            self.typ      = int(self.Header.Rat.rattype)
            self.info     = string.rstrip(self.Header.Info.info)
        elif self.version == 1.0:
            print "Warning: Old RAT format (convert!?)"
            e='<'
            if self.xdrflag == 1: e='>'
            lun = open(self.filename,'rb')      
            self.ndim = numpy.fromfile(file=lun, dtype=e+"i4",count=1)
            self.dim  = numpy.fromfile(file=lun, dtype=e+"i4",count=self.ndim)  
            self.var  = numpy.fromfile(file=lun, dtype=e+"i4",count=1)
            self.typ  = numpy.fromfile(file=lun, dtype=e+"i4",count=1)
            foo  = numpy.fromfile(file=lun, dtype=e+"i4",count=1)
            foo  = numpy.fromfile(file=lun, dtype=e+"i4",count=1)
            foo  = numpy.fromfile(file=lun, dtype=e+"i4",count=1)
            if self.xdrflag == 1: foo  = numpy.fromfile(file=lun, dtype=e+"i4",count=1)
            self.info  = string.rstrip(numpy.fromfile(file=lun, dtype="B",count=80).tostring())
            lun.close()
        else:
            print "ERROR: RAT version not supported" 
            raise IOError
        
        if   self.var == 1:  self.dtype = 'uint8'
        elif self.var == 2:  self.dtype = 'int16'
        elif self.var == 3:  self.dtype = 'int32'
        elif self.var == 4:  self.dtype = 'float32'
        elif self.var == 5:  self.dtype = 'float64'
        elif self.var == 6:  self.dtype = 'complex64'
        elif self.var == 9:  self.dtype = 'complex128'
        elif self.var == 12: self.dtype = 'uint16'
        elif self.var == 13: self.dtype = 'uint32'
        elif self.var == 14: self.dtype = 'int64'
        elif self.var == 15: self.dtype = 'uint64'
        else:
            print "ERROR: Array type not recognised!"
            raise IOError

            #--------------------------------------------------------------------------   
    def help(self):
        print
        print "FILE     : ",os.path.abspath(self.filename)
        if self.exists == True:
            print "SIZE     : ",self.dim[0:self.ndim]
            print "TYPE     : ",self.dtype
            print "INFO     : ",self.info
        else:
            print "--- empty file ---"
          
#--------------------------------------------------------------------------   
    def version(self):
        Header = RatHeader()
        magiclong  = Header.Rat.magiclong
        magicreal  = 0

        lun = open(self.filename,'rb')
        magicreal = numpy.fromfile(file=lun, dtype="i4",count=1)
        lun.close()
        
        if magicreal != magiclong:    # Check if maybe we have a RAT V1 File...
            lun = open(self.filename,'rb')
            ndim = numpy.fromfile(file=lun, dtype="<i4",count=1)
            lun.close()
            xdrflag = 0
            if ndim < 0 or ndim > 9:
                ndim = ndim.byteswap()  
                xdrflag = 1
            if ndim < 0 or ndim > 9:
                print "ERROR: format not recognised!"
                return False, False  
            version = 1.0
        else:                                                    #-------------- Yeah, RAT 2.0 found
            lun = open(self.filename,'rb')
            foo = numpy.fromfile(file=lun, dtype="int32",count=1)
            foo = numpy.fromfile(file=lun, dtype="float32",count=1)
            lun.close()
            version = foo[0]
            xdrflag = 0
            
        return version, xdrflag
