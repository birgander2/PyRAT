from __future__ import print_function
import ctypes
import os
import string
import numpy as np
import xml.etree.ElementTree as ET
from lxml import etree as xml_tools
import argparse
import glob
import logging
import shutil



def rrat(filename, **kwargs):
    """
    Reads an entire RAT file, returns it as a np ndarray variable
    """
    File = RatFile(filename)
    array = File.read(**kwargs)
    del File
    return array


def srat(filename, array, **kwargs):
    """Write a numpy ndarray into a RAT file."""
    File = RatFile(filename)
    File.write(array, **kwargs)
    del File


class RatHeaderRat(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("magiclong", ctypes.c_int),
                ("version", ctypes.c_float),
                ("ndim", ctypes.c_int),
                ("nchannel", ctypes.c_int),
                ("dim", ctypes.c_int * 8),
                ("var", ctypes.c_int),
                ("sub", ctypes.c_int * 2),
                ("rattype", ctypes.c_int ),
                ("reserved", ctypes.c_int * 9)]

    def __init__(self):
        self.magiclong = 844382546
        self.version = 2.0
        self.sub = (ctypes.c_int * 2)(1, 1)


class RatHeaderInfo(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("info", ctypes.c_char * 100)]


class RatHeaderGeo(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("projection", ctypes.c_short),
                ("ps_east", ctypes.c_double),
                ("ps_north", ctypes.c_double),
                ("min_east", ctypes.c_double),
                ("min_north", ctypes.c_double),
                ("zone", ctypes.c_short),
                ("hemisphere", ctypes.c_short),
                ("long0scl", ctypes.c_double),
                ("max_axis_ell", ctypes.c_double),
                ("min_axis_ell", ctypes.c_double),
                ("dshift_tx", ctypes.c_double),
                ("dshift_ty", ctypes.c_double),
                ("dshift_tz", ctypes.c_double),
                ("dshift_rx", ctypes.c_double),
                ("dshift_ry", ctypes.c_double),
                ("dshift_rz", ctypes.c_double),
                ("dshift_scl", ctypes.c_double),
                ("dshift_info", ctypes.c_char * 64),
                ("reserved", ctypes.c_byte * 18)]


class RatHeaderEmpty(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("reserved", ctypes.c_int * 25)]


class RatHeader(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("Rat", RatHeaderRat),
                ("Info", RatHeaderInfo),
                ("Geo", RatHeaderGeo),
                ("Stat", RatHeaderEmpty),
                ("Reserved1", RatHeaderEmpty),
                ("Reserved2", RatHeaderEmpty),
                ("Reserved3", RatHeaderEmpty),
                ("Reserved4", RatHeaderEmpty),
                ("Reserved5", RatHeaderEmpty)]

    def __init__(self):
        self.Rat = RatHeaderRat()


class RatFile():
    """
    Class for manipulating RAT formatted files. Work in progress, many features still missing...
    """

    def __init__(self, filename):
        self.filename = filename
        self.Header = RatHeader()
        if os.path.exists(self.filename):
            self.version, self.xdrflag = self.version()
            self.read_header()
            self.exists = True
        else:
            self.exists = False

    def read(self, **kwargs):
        """
        block = (rgstart, azstart, rgsize, azsize)
        """
        if self.version == 2.0:
            lun = open(self.filename, 'rb')
            if 'block' in kwargs:
                block = kwargs['block']
                if self.pixel_il:
                    blockshape = self.dim[::-1]+[1]*(2-self.ndim)
                    blockshape[0] = block[3]
                    arr = np.ndarray(shape=blockshape, dtype=self.dtype)
                    lun.seek(1000 + arr.itemsize * block[1] * blockshape[1] * self.nchannel)
                    lun.readinto(arr.data)
                    arr = arr[:, block[0]:block[0] + block[2], ...]
                else:
                    blockshape = self.dim[::-1]
                    blockshape[-2] = block[3]
                    ny = self.dim[1]
                    nx = self.dim[0]
                    band = np.ndarray(shape=blockshape[-2:], dtype=self.dtype)
                    blockshape[-1] = block[2]
                    arr = np.ndarray(shape=blockshape, dtype=self.dtype)
                    arr = arr.reshape((self.nchannel, block[3], block[2]))
                    seek_ch = nx * ny * arr.itemsize
                    for ch in range(self.nchannel):
                        lun.seek(1000 + ch * seek_ch + arr.itemsize * block[1] * nx)
                        lun.readinto(band)
                        arr[ch, ...] = band[:, block[0]:block[0] + block[2]]
                    arr = np.reshape(arr, blockshape)
            else:
                lun = open(self.filename, 'rb')
                lun.seek(1000)
                arr = np.ndarray(shape=self.dim[0:self.ndim][::-1], dtype=self.dtype)
                lun.readinto(arr.data)
            lun.close()
            return arr
        elif self.version == 1.0:
            lun = open(self.filename, 'rb')
            if 'block' in kwargs:
                block = kwargs['block']
                if self.pixel_il:
                    blockshape = self.dim[::-1]
                    blockshape[0] = block[3]
                    arr = np.ndarray(shape=blockshape, dtype=self.dtype)
                    lun.seek(104 + 4 * self.ndim + 4 * self.xdrflag + arr.itemsize * block[1] * blockshape[1] * self.nchannel)
                    lun.readinto(arr.data)
                    arr = arr[:, block[0]:block[0] + block[2], ...]
                else:
                    print("ERROR: block selection not supported for band interleaved images")
                    print("       in combination with RAT format version 1...")
                    arr = None
            else:
                lun.seek(104 + 4 * self.ndim + 4 * self.xdrflag)
                arr = np.ndarray(shape=self.dim[::-1], dtype=self.dtype)
                lun.readinto(arr.data)
            if self.xdrflag == 1:
                arr = arr.byteswap()
            lun.close()
            return arr
        else:
            print("ERROR: RAT version not supported")
            raise IOError

            # --------------------------------------------------------------------------

    def write(self, arr, type=0, **kwargs):
        if 'header' in kwargs:
            self.Header = kwargs['header']

        self.dtype = arr.dtype
        self.var = self.dtype2var(self.dtype)
        self.dim = arr.shape[::-1]
        self.ndim = arr.ndim
        self.deshape()
        self.Header.Rat.ndim = self.ndim
        self.Header.Rat.dim = self.dim
        self.Header.Rat.var = self.var
        self.Header.Rat.rattype = type
        self.Header.Rat.nchannel = self.nchannel

        self.info = self.Header.Info.info.decode().rstrip()

        lun = open(self.filename, 'wb')
        lun.write(self.Header)
        # lun.write(arr.tostring())
        arr.tofile(lun)
        lun.close()
        #self.help()

    def write_header(self):
        lun = open(self.filename, 'wb')
        lun.write(self.Header)
        return lun

    def dtype2var(self, dtype):
        if dtype == 'uint8' or dtype == 'bool':
            var = 1
        elif dtype == 'int16':
            var = 2
        elif dtype == 'int32':
            var = 3
        elif dtype == 'float32':
            var = 4
        elif dtype == 'float64':
            var = 5
        elif dtype == 'complex64':
            var = 6
        elif dtype == 'complex128':
            var = 9
        elif dtype == 'uint16':
            var = 12
        elif dtype == 'uint32':
            var = 13
        elif dtype == 'int64':
            var = 14
        elif dtype == 'uint64':
            var = 15
        else:
            var = 0
        return var

    def read_header(self):
        if self.version == 2.0:
            lun = open(self.filename, 'rb')
            lun.readinto(self.Header)
            lun.close()
            self.dim = [int(foo) for foo in self.Header.Rat.dim]
            self.ndim = int(self.Header.Rat.ndim)
            self.nchannel = int(self.Header.Rat.nchannel)
            self.deshape()
            self.var = int(self.Header.Rat.var)
            self.sub = [int(foo) for foo in self.Header.Rat.sub]
            self.typ = int(self.Header.Rat.rattype)
            #self.info = string.rstrip(self.Header.Info.info)
            self.info = self.Header.Info.info.rstrip()
        elif self.version == 1.0:
            # print("Warning: Old RAT format (convert!?)")
            e = '<'
            if self.xdrflag == 1: e = '>'
            lun = open(self.filename, 'rb')
            self.ndim = int(np.fromfile(file=lun, dtype=e + "i4", count=1))
            self.dim = np.fromfile(file=lun, dtype=e + "i4", count=self.ndim).tolist()
            self.deshape()
            # self.nchannel = int(np.prod(self.dim[2:self.ndim]))
            self.var = int(np.fromfile(file=lun, dtype=e + "i4", count=1))
            self.typ = int(np.fromfile(file=lun, dtype=e + "i4", count=1))
            foo = int(np.fromfile(file=lun, dtype=e + "i4", count=1))
            foo = int(np.fromfile(file=lun, dtype=e + "i4", count=1))
            foo = int(np.fromfile(file=lun, dtype=e + "i4", count=1))
            if self.xdrflag == 1: foo = np.fromfile(file=lun, dtype=e + "i4", count=1)
            self.info = np.fromfile(file=lun, dtype="B", count=80).tostring().rstrip().decode()
            lun.close()
        else:
            raise IOError("RAT version not supported")

        if self.var == 1:
            self.dtype = 'uint8'
        elif self.var == 2:
            self.dtype = 'int16'
        elif self.var == 3:
            self.dtype = 'int32'
        elif self.var == 4:
            self.dtype = 'float32'
        elif self.var == 5:
            self.dtype = 'float64'
        elif self.var == 6:
            self.dtype = 'complex64'
        elif self.var == 9:
            self.dtype = 'complex128'
        elif self.var == 12:
            self.dtype = 'uint16'
        elif self.var == 13:
            self.dtype = 'uint32'
        elif self.var == 14:
            self.dtype = 'int64'
        elif self.var == 15:
            self.dtype = 'uint64'
        else:
            print("ERROR: Array type not recognised!")
            raise IOError

#--------------------------------------------------------------------------
    def deshape(self):
        self.dim = self.dim[0:self.ndim]
        p1 = int(np.prod(self.dim[-2:]))
        p2 = int(np.prod(self.dim[0:-2]))
        p3 = int(np.prod(self.dim[0:2]))
        p4 = int(np.prod(self.dim[2:]))
        self.nchannel = min([p1, p2, p3, p4])
        idx = np.argmin([p1, p2, p3, p4])
        if idx == 0 or idx == 3:
            self.pixel_il = False
        else:
            self.pixel_il = True

    def help(self):
        print()
        print("FILE     : ", os.path.abspath(self.filename))
        if self.exists == True:
            print("SIZE     : ", self.dim[0:self.ndim])
            print("TYPE     : ", self.dtype)
            print("INFO     : ", self.info)
        else:
            print("--- empty file ---")

            #--------------------------------------------------------------------------

    def version(self):
        Header = RatHeader()
        magiclong = Header.Rat.magiclong
        magicreal = 0

        lun = open(self.filename, 'rb')
        magicreal = np.fromfile(file=lun, dtype="i4", count=1)
        lun.close()

        if magicreal != magiclong:  # Check if maybe we have a RAT V1 File...
            lun = open(self.filename, 'rb')
            ndim = np.fromfile(file=lun, dtype="<i4", count=1)
            lun.close()
            xdrflag = 0
            if ndim < 0 or ndim > 9:
                ndim = ndim.byteswap()
                xdrflag = 1
            if ndim < 0 or ndim > 9:
                raise IOError("format not recognised!")
                return False, False
            version = 1.0
        else:  #-------------- Yeah, RAT 2.0 found
            lun = open(self.filename, 'rb')
            foo = np.fromfile(file=lun, dtype="int32", count=1)
            foo = np.fromfile(file=lun, dtype="float32", count=1)
            lun.close()
            version = foo[0]
            xdrflag = 0

        return version, xdrflag


class Xml2Py(object):
    """Turns a STEP XML document (e.g. processing parameters) into a structure.

    Example usage:
    pp = Xml2Py(<path to pp file>)
    pp.v0
    (Python prints 89.1)
    pp.v0?
    (Type is numpy double)
    pp.r?
    (Numpy ndarray of doubles)

    Currently does not support complex data, pointer arrays and structure
    arrays.

     :author: Marc Jaeger

     :param fileName: Path to XML document.
     :type fileName: string

     :param elem: For internal use only! Do not use.
     :type elem: XML element object

    """

    def __init__(self, fileName, elem=None):
        if (elem == None):
            tree = ET.parse(fileName)
            root = tree.getroot()
            elem = root[0]

        self.xml2struct(elem)


    def xml2struct(self, elem):
        for f in elem:
            self.xml2field(f)


    @staticmethod
    def xml2val(elem):
        typElem = elem.find('datatype')
        dType = typElem.text
        dDim = typElem.attrib['length']
        dDim = np.asarray([np.long(d) for d in dDim.split()])[::-1]
        dLen = np.prod(dDim)

        valElem = elem.find('value')
        if (dType == 'pointer'):
            return Xml2Py.xml2val(valElem.find('parameter'))

        if (dType == 'struct'):
            return Xml2Py(None, valElem[0])

        val = elem.find('value').text
        if (dLen > 1):
            val = val.strip('[]').split(',')

        conv = {'int': np.int, 'long': np.long, 'float': np.float, 'double': np.double, 'string': lambda s:s}
        try:
            if (dLen > 1):
                val = np.asarray([conv[dType](v) for v in val]).reshape(dDim)
            else:
                val = conv[dType](val)
        except KeyError:
            print('WARNING: Unsupported data type {} in field {}! Ignoring...'.format(dType, 'n.a.'))
            return None

        return val



    def xml2field(self, elem):
        name = elem.attrib['name']
        setattr(self,name,self.xml2val(elem))




class Py2Xml(object):
    def __init__(self, template):
        self.__dict__['tree'] = ET.parse(template)
        self.__dict__['attrs'] = {}

        obj = self.tree.getroot().find('object')
        if (obj == None):
            raise ValueError('Expected an "object" element below the root!')

        self.get_attrs(obj)


    def get_attrs(self, obj, prefix=[]):
        for f in obj.iter('parameter'):
            name = prefix+[f.attrib['name']]

            while (f.find('./value/parameter') is not None):
                f = f.find('./value/parameter')

            v = f.find('value')
            t = f.find('datatype')

            if (v.find('object') is not None):
                self.get_attrs(v.find('object'),name)
                continue

            self.attrs['_'.join(name)] = (v,t,f)


    def __getattr__(self,name):
        if (name in self.__dict__):
            return self.__dict__[name]
        v,t,f = self.__dict__['attrs'][name]
        return Xml2Py.xml2val(f)


    def __setattr__(self,name,value):
        if (name in self.__dict__):
            self.__dict__[name] = value
        v,t,f = self.attrs[name]
        return self.val2xml(v,t,value)

#,long: (str,'long')
    @staticmethod
    def val2xml(v, t, value):
        cdict = {str: (str, 'string'),float:(str,'double'),int: (str,'long'),complex:(lambda z:'({},{})'.format(z.real,z.imag),'complex')}

        cdict[np.int32] = cdict[int]
        cdict[np.int64] = cdict[int]
        cdict[np.float32] = (str,'float')
        cdict[np.float64] = cdict[float]
        cdict[np.complex64] = cdict[complex]

        # arrays with this shape in xml, e.g. [x1,y1,...,xN,yN] :
        # <datatype length="2 4">double</datatype>
        # have this shape (4, 2) in python/numpy
        # therefore value.shape is reversed : value.shape[::-1]
        if (isinstance(value, np.ndarray)):
            t.attrib['length'] = ' '.join([str(l) for l in value.shape[::-1]])
            vtype = type(value.flat[0])
            t.text = cdict[vtype][1]
            v.text = '['+', '.join([cdict[vtype][0](val) for val in value.flat])+']'
        else:
            try:
                if isinstance(value,str):
                    raise ValueError()
                t.attrib['length'] = str(len(value))
                t.text = cdict[type(value[0])][1]
                v.text = '['+', '.join([cdict[type(value[0])][0](val) for val in value])+']'
            except:
                t.attrib['length'] = '1'
                t.text = cdict[type(value)][1]
                v.text = cdict[type(value)][0](value)


    def set_attr(self, name, value, prefix=[]):

        if isinstance(value,dict):
            for k in value:
                self.set_attr(k,value[k],prefix+[name])
            return

        try:
            v,t,f = self.attrs['_'.join(prefix+[name])]
            self.val2xml(v, t, value)
        except KeyError:
            pass


    def update(self, obj):
        d = obj
        if (hasattr(obj,'__dict__')):
            d = obj.__dict__

        if not isinstance(d,dict):
            raise ValueError('Expected a dictionary or an object with a __dict__ attribute!')

        for k in d:
            self.set_attr(k,d[k])

    def write(self, filename, overwrite=False):
        if not overwrite:
            if os.path.isfile(filename):
                raise FileExistsError
        self.tree.write(filename)

    def __str__(self):
        return str(self.__dict__)
