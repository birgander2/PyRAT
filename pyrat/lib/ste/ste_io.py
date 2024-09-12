""""
Input-output Module
===================

The following module provides the possibility to work with ``RAT`` files, widely
used by DLR-HR institute and some others.

:author: Andreas Reigber <andreas.reigber@dlr.de>
:author: Anton Heister <anton.heister@dlr.de>

"""

import ctypes
import os
import copy

import numpy as np
import mmap
from functools import reduce

import warnings

try:
    import lxml.etree as ET
except ImportError:
    warnings.warn('Failed to import ET from "lxml.etree". The Py2XML class will not work.', ImportWarning)

try:
    from pkg_resources import resource_string
    from mako.template import Template
except ImportError:
    warnings.warn('Failed ot import pkg_resources and/or mako. Writing geo-Envi-headers for RAT files will not work.',
                  ImportWarning)

red = "\033[91m"
endc = "\033[0m"


def rrat(filename, **kwargs):
    """Read an entire RAT file, return it as a numpy array."""
    rat_file = RatFile(filename)
    return rat_file.read(**kwargs)

def mrrat(filename, **kwargs):
    """Read an entire RAT file, return it as a memory map to the numpy array.

    Convinient for reading big files especially for blockwise processing. The
     function works faster than rrat, but the disadvantage is that the array is
     read-only.
    """
    rat_file = RatFile(filename)
    return rat_file.mread(**kwargs)

def srat(filename, array, **kwargs):
    """Write a numpy ndarray into a RAT file."""
    rat_file = RatFile(filename)
    rat_file.write(array, **kwargs)


class RatHeaderRat(ctypes.Structure):
    """
    A base class to store RAT header's attributes.

    The class is a child of ``ctypes.Structure``. This is made to write and read
    RAT headers correctly and easily preserving the size of each field in bytes.

    :param magiclong: used for RAT version control
    :type magiclong: int
    :param version: RAT version (currently 1.0 and 2.0 versions are available)
    :type version: float
    :param ndim: number of array's dimension
    :type ndim: int
    :param nchannel: number of channels
    :type nchannel: int
    :param shape: the shape of the data array
    :type shape: list
    :param var: specifies data type according to IDL's convention. More
      information here <http://www.exelisvis.com/docs/IDL_Data_Types.html>
    :type var: int
    :param sub: to be implemented
    :type sub: int
    :param rattype: to be implemented
    :type rattype: int

    .. note:: Due to some degree of redundancy in RatHeader  (there are 3
      connected attributes: ``shape``, ``ndim``, ``nchannel``), the priority in
      determining the value of ``ndim`` and ``nchannel`` parameters is given to
      ``shape`` parameter. That means that ``ndim`` is calculated as
      ``len(shape)`` and ``nchannel`` can be given as ``**kwargs``, if not given
      then it is equal either to ``product(shape[2::])``  if ``ndim > 2`` or to ``1`` if
      ``ndim < 2``.
    """
    _pack_ = 1
    _fields_ = [("magiclong", ctypes.c_int),
                ("version", ctypes.c_float),
                ("ndim", ctypes.c_int),
                ("nchannel", ctypes.c_int),
                ("idl_shape", ctypes.c_int * 8),
                ("var", ctypes.c_int),
                ("sub", ctypes.c_int * 2),
                ("rattype", ctypes.c_int),
                ("reserved", ctypes.c_int * 9)]


    def __init__(self, **kwargs):
        """Create an object of RatHeaderRat class.

        Create RatHeaaderRat instance. Due to redundancy of ``RAT`` header
        (there are 3 connected attributes: ``shape``, ``ndim``, ``nchannel``).
        Until the latter 2 are not specified explicitly, their value is
        calculated based on ``shaped``.

        To specify the data type of a ``RAT`` file one may use either ``dtype``
        keyword with numpy data type or a string as a value, or a ``var``
        keyword according to IDL's regulations.

        **Keywords**:
          :param array: a numpy array which is used to define/override ``shape`` and ``dtype``.
          :type shape: list
          :param shape: a shape of a ``RAT`` file.
          :type shape: list
          :param ndim: number of dimensions, if not given, then parsed from
            ``shape``
          :type ndim: int
          :param nchannel: number of channels, if not given, then parsed from
            ``shape``.
          :type nchannel: int
          :param var: a data type according to IDL's regulation
          :type var: int
          :param dtype: a data type of a ``RAT`` file.
          :type dtype: numpy.dtype or a string
          :param dtype: a data type of a ``RAT`` file.
          :type dtype: numpy.dtype or a string
          :param sub: ...
          :param sub: list
          :param rattype: type of a ``RAT`` file
          :type rattype: int

        """

        if 'array' in kwargs:
            kwargs['shape'] = kwargs['array'].shape
            kwargs['dtype'] = kwargs['array'].dtype

        self.magiclong = 844382546
        self.version = 2.0
        # checking kwargs
        if 'shape' in kwargs:
            # reverse the shape so it corresponds to RAT convention
            self.idl_shape = (ctypes.c_int * 8)(*kwargs['shape'][::-1])
            if 'ndim' in kwargs:
                self.ndim = ctypes.c_int(kwargs['ndim'])
            else:
                self.ndim = ctypes.c_int(len(kwargs['shape']))

            if 'nchannel' in kwargs:
                self.nchannel = ctypes.c_int(kwargs['nchannel'])
            else:
                if len(kwargs['shape']) <= 2:
                #    self.nchannel = ctypes.c_int(0)
                    self.nchannel = ctypes.c_int(1)
                else:
                #    self.nchannel = ctypes.c_int(kwargs['shape'][-1])
                    self.nchannel = ctypes.c_int(np.prod(kwargs['shape'][2::]))
        if 'var' in kwargs:
            self.var = ctypes.c_int(kwargs['var'])
        if 'rattype' in kwargs:
            self.rattype = ctypes.c_int(kwargs['rattype'])
        if 'sub' in kwargs:
            self.sub = (ctypes.c_int * 2)(*kwargs['sub'])
        else:
            self.sub = (ctypes.c_int * 2)(1, 1)
        if 'dtype' in kwargs:
            if type(kwargs['dtype']) == type:
                data_type = np.dtype(kwargs['dtype']).name
            elif type(kwargs['dtype'] == str):
                data_type = kwargs['dtype']
            self.var = get_var(data_type)

class RatHeaderInfo(ctypes.Structure):
    """Contains a 100 character line for RatFile description.
    """
    _pack_ = 1
    _fields_ = [("info", ctypes.c_char * 100)]


class RatHeaderGeo(ctypes.Structure):
    """Contains positioning information.

    :param projection:
    :type projection:
    :param ps_east:
    :type ps_east:
    :param ps_north:
    :type ps_north:
    :param min_east:
    :type min_east:
    :param min_north:
    :type min_north:
    :param zone:
    :type zone:
    :param hemisphere:
    :type hemisphere:
    :param long0scl:
    :type long0scl:
    :param max_axis_ell:
    :type max_axis_ell:
    :param min_axis_ell:
    :type min_axis_ell:
    :param dshift_tx:
    :type dshift_tx:
    :param dshift_ty:
    :type dshift_ty:
    :param dshift_tz:
    :type dshift_tz:
    :param dshift_rx:
    :type dshift_rx:
    :param dshift_ry:
    :type dshift_ry:
    :param dshift_rz:
    :type dshift_rz:
    :param dshift_scl:
    :type dshift_scl:
    :param dshift_info:
    :type dshift_info:
    """
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
    """Contains RatHeaderRat, RatHeaderInfo and RatHeaderGeo.

        Create RatHeaaderRat instance. Due to redundancy of ``RAT`` header
        (there are 3 connected attributes: ``shape``, ``ndim``, ``nchannel``).
        Until the latter 2 are not specified explicitly, their value is
        calculated based on ``shaped``.

        To specify the data type of a ``RAT`` file one may use either ``dtype``
        keyword with numpy data type or a string as a value, or a ``var``
        keyword according to IDL's regulations.

        **Keywords**:
          :param shape: a shape of a ``RAT`` file.
          :type shape: list
          :param ndim: number of dimensions, if not given, then parsed from
            ``shape``
          :type ndim: int
          :param nchannel: number of channels, if not given, then parsed from
            ``shape``.
          :type nchannel: int
          :param var: a data type according to IDL's regulation
          :type var: int
          :param dtype: a data type of a ``RAT`` file.
          :type dtype: numpy.dtype or a string
          :param sub: ...
          :param sub: list
          :param rattype: type of a ``RAT`` file
          :type rattype: int

        """
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

    def __init__(self, **kwargs):
        self.Rat = RatHeaderRat(**kwargs)


class RatFile():
    """    Class for manipulating RAT formatted files."""

    def __init__(self, filename):
        """Initialize a RAT file.

        If the file exists reads the file's header, if not then initializes
        the empty ``RatFile`` **instance** with a given filename.

        :param filename: either a RAT filename (with *.rat extension) if the
          file is in current working directory or an absolute path of the file.
        :type filename: string
        :return: RatFile instance

        .. note:: If the file doesn't exist the function creates a new
          ``RatFile`` instance with an empty header, not an empty file.
        """
        self.filename = filename
        self.Header = RatHeader()
        # the shape of numpy array
        self.shape = ()
        try:
            self.version, self.xdrflag = self.get_version()
            self.read_header()
            self.shape = self._get_shape()
            self.ndim = len(self.shape)
            self.dtype = self._get_dtype()
            self.info = self.Header.Info.info.decode()
            self.var = int(self.Header.Rat.var)
            self.nchannel = int(self.Header.Rat.nchannel)
            self.exists = True
        except (IOError, IndexError):
            self.exists = False


    @classmethod
    def _ioerror(cls, msg):
        print(red + msg + endc)
        raise IOError(msg)

    def create(self, shape=None, header=None, **kwargs):
        """Create an empty ``RAT`` file and write a RAT header into it.

        Create an empty ``rat`` file with given parameters and write a RAT
        header. Either ``RatHeader`` instance or ``shape`` list should be
        given as an argument. If both are given, then the value of ``shape``
        overrides the ``RatHeaderRat.shape`` value. It is necessary to specify
        a ``dtype`` keyword if the proper datatype wasn't specified in
        header's ``RatHeaderRat.var`` attribute).

        :param shape: the shape of the data to store in *.rat file
        :type shape: list
        :param header: a rat header
        :type header: RatHeaderRat
        :keyword dtype: data type
        :type dtype: type or string

        :return: None

        :raises: IOError

        .. note::

          1) The ``dtype`` keyword maybe given either as a string
             (i.e. 'int16') or as a numpy ``type`` instance (i.e. ``np.int16``).

          2) The function creates an empty binary file which is a sparse file
             for Linux and Windows, OS X however doesn't support sparse files
             so it will be just an empty file there.

        """
        # raise an error if neither header nor shape is given as an arg
        if (header is None) and (shape is None):
            self._ioerror('Please, specify either shape or header!')

        if header is not None:
            self.Header = header

        if shape is not None:
            self.Header.Rat.idl_shape = (ctypes.c_int * 8)(*shape[::-1])

        if 'dtype' in kwargs:
            if type(kwargs['dtype']) == type:
                data_type = np.dtype(kwargs['dtype']).name
            elif type(kwargs['dtype'] == str):
                data_type = kwargs['dtype']
            self.Header.Rat.var = get_var(data_type)

        self.shape = self._get_shape()
        self.dtype = self._get_dtype()
        self.Header.Rat.ndim = ctypes.c_int(len(self.shape))
        self.Header.Rat.nchannel = ctypes.c_int(int(np.prod(self.shape[2:])))

        if 'rattype' in kwargs:
            self.Header.Rat.rattype = ctypes.c_int(kwargs['rattype'])

        # calculate the needed size of an empty file
        n_bytes = reduce(lambda x, y: x * y, self.shape) * self.dtype.itemsize

        # write the Header and truncate the file
        with open(self.filename, 'wb') as lun:
            lun.write(self.Header)
            #lun.truncate(n_bytes)
            self.exists = True

        return self


    def write(self, arr=[], **kwargs):
        """Write either a whole data array or a block of data into a rat file.

        The following usages are available:

        Write a whole data array either with or without the ``header`` keyword.
        In latter case the parameters of the array (shape and dtype), are parsed
        into an empy ``header`` instance and written to the file. In this case
        it's also possible to explicitly specify ``rattype`` keyword.

        When specifing the ``header`` keyword it is also possible to specify
        ``shape``, ``dtype`` and ``rattype`` keywords so the values of these
        parameters contained in ``header`` will be overwritten.

        When using a ``block`` keyword it is supposed that the header was
        already written (the entire shape of the ``RAT`` file should be known
        prior block writting).

        To write only a ``header`` one should use ``rat.write(header=header)``
        command, where ``rat`` is a ``RatFile`` instance and ``header`` is a
        ``RatHeader`` instance.

        Using a ``block`` keyword for 2D arrays the offset in both axes is
        possible (i.e. if the shape specified in header is [100,200], then a
        block can be equal to [20,50, 130, 200]). In other cases only an offset
        along the first axis is supported. The data type of the ``arr`` should
        correspond to the one given in previously written ``header``, otherwise
        an error will be raised.

        :param arr: array to be stored in rat file
        :type arr: numpy.ndarray

        **Keywords**:
          :param block: a position, where to write an array; has the following
            form ``[start_1, stop_1, ..., ..., start_N, stop_N]``, where ``N``
            is a number of array's dimensions.
          :type block: list
          :param header: RAT header
          :type header: RatHeaderRat
          :param shape: shape to overwrite the shape of the header
          :type shape: list
          :param dtype: data type to be written into header
          :param rattype: numpy.dtype or a string
          :param rattype: specifies RAT file type
          :type rattype: int


        :raises: IOError

        """

        arr = np.asarray(arr)
        if arr.ndim == 0:
            arr = arr.reshape(1)

        # block writing
        if 'block' in kwargs:
            # check if datatype of arr and of 'var' are the same
            self._check_dtypes(arr)

            block = kwargs['block']
            # check if the block var meets the requirements
            block = self._check_block(block, arr=arr)

            if 'header' in kwargs:
                self._ioerror('The header should have been written prior to block '
                              'writting!')

            # for 2D arrays an offset for both axes is allowed
            if arr.ndim == 2:
                # create zero-filled numpy array of self.shape
                # and fill it with arr.data corresponding to block array values
                temp = np.zeros(
                    [arr.shape[0], self.shape[1]], dtype=arr.dtype)
                temp[:, block[2]:block[3]] = arr
                arr = temp
                offset = 1000 + self.shape[1] * block[0] * arr.itemsize

            else:
                offset = 1000 + np.ravel_multi_index(
                    block[::2], self.shape) * arr.itemsize

            with open(self.filename, 'r+') as lun:
                # The header should be written prior block writting
                lun.seek(offset)
                arr.tofile(lun)
                lun.flush()
                self.exists = True
            return

        # no block writing
        else:
            if 'header' in kwargs:
                self.Header = kwargs['header']
                # modify existing header if needed
                if 'shape' in kwargs:
                    self.Header.Rat.idl_shape = (
                        ctypes.c_int * 8)(*kwargs['shape'][::-1])
                    self.shape = tuple(kwargs['shape'])
                if 'dtype' in kwargs:
                    if type(kwargs['dtype']) == type:
                        data_type = np.dtype(kwargs['dtype']).name
                    elif type(kwargs['dtype'] == str):
                        data_type = kwargs['dtype']
                    self.Header.Rat.var = get_var(data_type)
                if 'rattype' in kwargs:
                    self.Header.Rat.rattype = kwargs['rattype']
                # check if datatypes of array and header are equal
                if arr.size > 0:
                    self._check_dtypes(arr)

            elif arr.size > 0:
                # parse array parameters to the Header
                self.Header = RatHeader(
                    shape=(arr.shape if 'shape' not in kwargs else kwargs['shape']),
                    var=np.ctypeslib.array(get_var(arr.dtype)))
                if 'rattype' in kwargs:
                    self.Header.Rat.rattype = kwargs['rattype']

            else:
                self._ioerror('Specify a header, an array or both!')

            self.dtype = self._get_dtype()
            self.shape = self._get_shape()
            self.Header.Rat.ndim = ctypes.c_int(len(self.shape))
            self.Header.Rat.nchannel = ctypes.c_int(int(np.prod(self.shape[2:])))

            n_bytes_total = (
            1000 + reduce(lambda x, y: x * y, self.shape) * self.dtype.itemsize)

            with open(self.filename, 'wb') as lun:
                lun.write(self.Header)
                if arr.size > 0:
                    arr.tofile(lun)
                self.exists = True
                if lun.tell() > n_bytes_total:
                    warnings.warn("The size of the RAT file exceed! The array "
                                  "is written outside the header's dimensions!")
                lun.flush()
            return
    # --------------------------------------------------------------------------
    def read(self, **kwargs):
        """Read the data from ``RAT`` file as a numpy array.

        Works both with ``RAT`` 1.0 and 2.0 files, allows to read the data in
        blocks along all the axes.

        **Keywords**:
          :param block: the block of data to read;
            ``block=[start_1, stop_1, ..., ..., start_N, stop_N]``, where ``N``
            is a number of axes.

        :return: numpy.ndarray

        :raises: IOError if the file doesn't exist, if ``RAT`` version is not
          recognized and when ``block`` doesn't correspond to the shape of
          ``RAT`` header.
        """
        if self.exists == False:
            self._ioerror('ERROR: The file "%s" does not exist'%self.filename)
        if 'block' in kwargs:
            block = kwargs['block']
            # check if the block var meets the requirements
            block = self._check_block(block)
        else:
            block = np.zeros(2 * len(self.shape), dtype=np.int32)
            block[1::2] = self.shape

        ind = tuple(map(
            lambda x, y: slice(x, y, None), block[::2], block[1::2]))

        if self.version == 2.0:
            offset = 1000
        elif self.version == 1.0:
            offset = int(104 + 4 * self.Header.Rat.ndim + 4 * self.xdrflag)
        else:
            self._ioerror('ERROR: RAT version not supported')

        with open(self.filename, 'rb') as lun:
            mm = mmap.mmap(
                lun.fileno(), length=0, access=mmap.ACCESS_READ)

            arr = (np.ndarray.__new__(np.ndarray, self.shape, dtype=self.dtype,
                                      buffer=mm, offset=offset)[ind])
            if self.xdrflag == 1:
                arr = arr.byteswap()
        arr_new = np.zeros(shape=arr.shape, dtype=arr.dtype)
        arr_new[:] = arr
        return arr_new

    # --------------------------------------------------------------------------
    def mread(self, **kwargs):
        """Read the data from ``RAT`` file as a numpy array fast using a memory
        map (attention: the array is opened in read-only mode).

        Works both with ``RAT`` 1.0 and 2.0 files, allows to read the data in
        blocks along all the axes.

        **Keywords**:
          :param block: the block of data to read;
            ``block=[start_1, stop_1, ..., ..., start_N, stop_N]``, where ``N``
            is a number of axes.

        :return: numpy.ndarray

        :raises: IOError if the file doesn't exist, if ``RAT`` version is not
          recognized and when ``block`` doesn't correspond to the shape of
          ``RAT`` header.
        """
        if self.exists == False:
            self._ioerror('ERROR: The file is not found')
        if 'block' in kwargs:
            block = kwargs['block']
            # check if the block var meets the requirements
            block = self._check_block(block)
        else:
            block = np.zeros(2 * len(self.shape), dtype=np.int32)
            block[1::2] = self.shape

        ind = tuple(map(
            lambda x, y: slice(x, y, None), block[::2], block[1::2]))

        if self.version == 2.0:
            offset = 1000
        elif self.version == 1.0:
            offset = int(104 + 4 * self.Header.Rat.ndim + 4 * self.xdrflag)
        else:
            self._ioerror('ERROR: RAT version not supported')

        with open(self.filename, 'rb') as lun:
            mm = mmap.mmap(
                lun.fileno(), length=0, access=mmap.ACCESS_READ)

            arr = (np.ndarray.__new__(np.ndarray, self.shape, dtype=self.dtype,
                                      buffer=mm, offset=offset)[ind])
            if self.xdrflag == 1:
                arr = arr.byteswap()
        return arr

    def append(self, arr):
        """Append the ``RAT`` file with a given array along the first axis.

        :param arr: the array to write
        :type arr: numpy.ndarray

        :raises: warning if the dimensions specified in file's header are exceed
        """
        # self.read_header()
        # check if datatype of arr and of 'var' are the same
        self._check_dtypes(arr)

        if arr.ndim == len(self.shape)-1:
            arr = arr[np.newaxis,...]

        # check if the array's and header's shape correspond to each other
        if len(self.shape) > 1 and (self.shape[1:] != arr.shape[1:]):
            self._ioerror('The shape specified in the header %s and the shape of '
                          'the array %s don\'t correspond to each other!' % (str(self.shape[1:]),str(arr.shape[1:])))

        n_bytes_total = 1000 + reduce(lambda x, y: x * y, self.shape) * arr.itemsize

        with open(self.filename, 'r+b') as lun:
            lun.seek(0,2)
            arr.tofile(lun)
            if lun.tell() > n_bytes_total:
                warnings.warn("The size of the RAT file exceed! The array is "
                              "written outside the header's dimensions!")

        return


    # --------------------------------------------------------------------------
    def read_header(self):
        # reading the version
        self.version = self.get_version()[0]
        """Read ``RAT`` header; supports both ``RAT`` 1.0 and 2.0 versions."""
        if self.version == 2.0:
            with open(self.filename, 'rb') as lun:
                lun.readinto(self.Header)

        elif self.version == 1.0:
            warnings.warn('Old RAT v1.0 format!')
            if self.xdrflag == 1:
                data_type = '>i4'
                offset = 4 * 4
            else:
                data_type = '<i4'
                offset = 3 * 4

            with open(self.filename, 'rb') as lun:
                ndim = int(np.fromfile(file=lun, dtype=data_type, count=1))
                shape = np.fromfile(file=lun, dtype=data_type, count=ndim).tolist()
                shape = shape[::-1]  
                var = int(np.fromfile(file=lun, dtype=data_type, count=1))
                rattype = int(np.fromfile(file=lun, dtype=data_type, count=1))
                lun.seek(offset, 1)
                info = np.fromfile(file=lun, dtype="B", count=80).tostring().rstrip()

            # initialize the header
            self.Header = RatHeader(shape=shape, ndim=ndim, var=var,
                                    rattype=rattype)
            self.Header.Info.info = info

        else:
            self._ioerror('ERROR: RAT version not supported')

        self.dtype = self._get_dtype()

    #--------------------------------------------------------------------------

    def write_envi_header(self, info='', sensorType='DLR F-SAR'):
        hdrFile = self.filename+'.hdr'
        with open(hdrFile,'w') as f:
            f.write('ENVI\n')
            if (len(info) > 0):
                f.write('description     = {%s}\n' % info)
            f.write('samples         = %i\n' % self.Header.Rat.idl_shape[0])
            f.write('lines           = %i\n' % self.Header.Rat.idl_shape[1])
            f.write('bands           = %i\n' % np.maximum(self.Header.Rat.nchannel,1))
            f.write('header offset   = 1000\n')
            f.write('file type       = ENVI Standard\n')
            f.write('data type       = %i\n' % self.Header.Rat.var)
            f.write('interleave      = bsq\n')
            f.write('byte order      = 0\n')
            if (len(sensorType) > 0):
                f.write('sensor type     = %s\n' % sensorType)

    def write_geo_envi_header(self):
            tmpl = Template(resource_string(__package__+'.templates', 'envihdr.tpl'))
            envi_hdr = tmpl.render(file=self.filename, hdr=self.Header)
            with open(self.filename+'.hdr','w') as f:
                f.write(envi_hdr)

    #--------------------------------------------------------------------------

    def help(self):
        """Print basic imformation about ``RAT`` file.

        Print ``shape``, ``data type`` and ``info`` from file's header.
        If the file doesn't exist, then print ``empy file``.
        """
        print()
        print("FILE  : ", os.path.abspath(self.filename))
        if self.exists == True:
            print("SHAPE : ", self.shape)
            print("TYPE  : ", self._get_dtype())
            print("INFO  : ", self.Header.Info.info)
        else:
            print("--- empty file ---")

    #--------------------------------------------------------------------------

    def get_version(self):
        """Get the version of ``RAT`` file: 1.0 or 2.0."""
        Header = RatHeader()
        magiclong = Header.Rat.magiclong
        magicreal = 0

        with open(self.filename, 'rb') as lun:
            magicreal = np.fromfile(file=lun, dtype="i4", count=1)
        if magicreal != magiclong:  # Check if maybe we have a RAT V1 File...
            with open(self.filename, 'rb') as lun:
                ndim = np.fromfile(file=lun, dtype="<i4", count=1)
            xdrflag = 0
            if ndim < 0 or ndim > 9:
                ndim = ndim.byteswap()
                xdrflag = 1
            if ndim < 0 or ndim > 9:
                print(red + "ERROR: format not recognised!" + endc)
                return False, False
            version = 1.0
        else:  #-------------- Yeah, RAT 2.0 found
            with open(self.filename, 'rb') as lun:
                lun.seek(4)
                version = np.fromfile(file=lun, dtype="float32", count=1)[0]
            xdrflag = 0

        return version, xdrflag

    def _check_block(self, block, **kwargs):
        if len(block) == 4:                        # only 2D block provided
            block = list(block)
            dimlist = list(self.shape)
            dimlist[dimlist.index(max(dimlist))] = 0
            dimlist[dimlist.index(max(dimlist))] = 0
            for k, dim in enumerate(dimlist):
                if dim != 0:
                    block.insert(k*2, dim)
                    block.insert(k*2, 0)

        stop_more_than_shape = reduce(lambda x, y: x or y, (
            map(lambda x, y: (x > y), block[1::2], self.shape)))
        if stop_more_than_shape:
            self._ioerror('Value of block exceeds the array shape!')

        block = np.asarray(block)
        if block.dtype.kind not in ('i','u'):
            self._ioerror('Block extent must be given by integers!')

        if np.min(block) < 0:
            self._ioerror('The items in block must be nonnegative!')

        if 'arr' in kwargs:
            if len(block) // 2 != kwargs['arr'].ndim:
                self._ioerror('The dimensions of block do not correspond to the '
                              'dimensions of array!')

            block_not_shape = reduce(lambda x, y: x or y,
                                     map(lambda x, y, z: (x - y) != z,
                                         block[1::2], block[::2],
                                         kwargs['arr'].shape))
            if block_not_shape:
                self._ioerror('Length of block components %s does not correspond to'
                              ' the shape of the array %s!'%(str(block),str(kwargs['arr'].shape)))
        return block


    def _check_dtypes(self, arr):
        """Check if dtypes of given array and the one in header are equal"""
        if self.Header.Rat.var != get_var(arr.dtype).value:
            self._ioerror('The data type of the array to be written and '
                          'the one specified in file\'s header '
                          'don\'t correspond to each other!')

    def _get_dtype(self):
        """Get data type give ```Header.Rat.var``."""
        try:
            return np.dtype(dtype_dict[self.Header.Rat.var])
        except KeyError:
            self._ioerror('The data type is either not specified or '
                          'not supported!')

    def _get_shape(self):
        """Get numpy array shape from ``Header.Rat.shape`` that is IDL style"""
        shape = np.ctypeslib.as_array(self.Header.Rat.idl_shape)
        shape = shape[shape != 0]
        return tuple(shape[::-1])


def check_ratformat(filename):
    with open(filename, 'rb') as lun:
        magiclong = lun.read(4)
        if magiclong == b'RAT2':
            return True
        else:
            return False


def get_var(dtype):
    """Get ``RatHeaderRat.var`` value given ``dtype``."""
    var = [key for (key, value) in dtype_dict.items() if
           value == dtype]
    if var == []:
        RatFile._ioerror('The data type is not supported!')
    else:
        return ctypes.c_int(var[0])

# data type dictionary to net RAT's and np's data formats
dtype_dict = {1: 'uint8',
              2: 'int16',
              3: 'int32',
              4: 'float32',
              5: 'float64',
              6: 'complex64',
              9: 'complex128',
              12: 'uint16',
              13: 'uint32',
              14: 'int64',
              15: 'uint64'}





class Py2Xml(object):
    """
    Class to read/write DLR-STE XML documents (primarily used by the STEP processor for processing parameters)

    Allows python object-like access to parameters.

    For example, for STEP processing parameters (XML files in RGI/RGI-RDP):

    pp = Py2Xml('path_to_pp.xml')
    print(pp.v0)     # v0 is automatically converted to floating point
    print(pp.r[:10]) # first 10 values of range vector, also floating point

    The parameter object also allows dictionary-like access to handle problematic parameter names
    (which clash with python keywords). For example:

    print(pp['lambda']) # pp.lambda would be a syntax error
    print(pp['pass'])   # same as above
    print(pp['v0'])     # dictionary style access works for other parameters, too!

    The class provides full read/write support. Parameter values are changed by standard assignment
    and structures can be saved using the write method:

    pp.v0 = 100
    pp.write('path_to_new_pp.xml')

    """

    def __init__(self, root):
        if isinstance(root, str):
            self.__dict__['__root__'] = ET.parse(root).find('object')
        else:
            self.__dict__['__root__'] = root

        if self.__root__ is None:
            raise ValueError('Expected an "object" element below the root element!')

        self.__dict__['__iterend__'] = False


    def __getstate__(self):
        return ET.tostring(self.__root__, encoding=str)

    def __setstate__(self, root):
        self.__dict__['__root__'] = ET.fromstring(root)
        self.__dict__['__iterend__'] = False

    def copy(self):
        return Py2Xml(copy.deepcopy(self.__dict__['__root__']))

    def __getparam__(self, name):
        p = [p for p in self.__root__.iterchildren('parameter') if p.attrib['name'] == name]
        if len(p) != 1:
            raise AttributeError('Expected a unique match parameter name "%s", got %i matches.' % (name, len(p)))

        return [p[0].find(tag) for tag in ('remark', 'datatype', 'value')]


    def params(self):
        return [p.attrib['name'] for p in self.__root__.iterchildren('parameter')]

    @staticmethod
    def xml2val(v, t):
        type = t.text
        shape = t.attrib['length']
        shape = np.asarray([np.uint64(d) for d in shape.split()])[::-1]
        size = np.prod(shape)

        if (type == 'pointer'):
            p = v.find('parameter')
            return Py2Xml.xml2val(*[p.find(t) for t in ('value', 'datatype')])

        if (type == 'struct'):
            obj_arr = [Py2Xml(obj) for obj in v.iterchildren('object')]
            return obj_arr[0] if size <= 1 else obj_arr

        conv = {'int': int, 'long': int, 'float': np.float32, 'double': np.float64, 'string': lambda s: s}
        try:
            if size > 1:
                val = np.asarray([conv[type](v) for v in v.text.strip('[]').split(',')]).reshape(shape)
            else:
                val = conv[type](v.text)
        except KeyError:
            print('PY2XML WARNING: Unsupported data type "%s" encountered. Skipping!' % (type))
            return None

        return val

    @staticmethod
    def val2xml(v, t, value):
        cdict = {str: (str, 'string'), \
                 int: (str, 'long'), \
                 float: (str, 'double'), \
                 complex: (lambda z: '({},{})'.format(z.real, z.imag), 'complex')}
        cdict[np.uint8] = cdict[int]
        cdict[np.int32] = cdict[int]
        cdict[np.uint32] = cdict[int]
        cdict[np.int64] = cdict[int]
        cdict[np.uint64] = cdict[int]
        cdict[np.float32] = (str, 'float')
        cdict[np.float64] = cdict[float]
        cdict[np.complex64] = cdict[complex]
        cdict[bool] = cdict[str]

        if (t.text == 'pointer'):
            p = v.find('parameter')
            return Py2Xml.val2xml(*([p.find(t) for t in ('value', 'datatype')] + [value]))

        try:
            vsize = 1 if isinstance(value, str) else len(value)
        except TypeError:
            vsize = 1

        t.attrib['length'] = str(vsize)
        v.clear()
        if vsize == 1 and not isinstance(value, Py2Xml):
            t.text = cdict[type(value)][1]
            v.text = cdict[type(value)][0](value)
        elif all([isinstance(v, Py2Xml) for v in value]):
            t.text = 'struct'
            for obj in value:
                v.append(copy.deepcopy(obj.__root__))
        else:
            if isinstance(value, np.ndarray):
                t.attrib['length'] = ' '.join([str(l) for l in value.shape[::-1]])
                value = value.flat
            vtype = type(value[0])
            t.text = cdict[vtype][1]
            v.text = '[' + ', '.join([cdict[vtype][0](val) for val in value]) + ']'

    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        if key == 0:
            return self
        r, t, v = self.__getparam__(key)
        return Py2Xml.xml2val(v, t)

    def __getitem__(self, key):
        return self.__getattr__(key)

    def __setattr__(self, name, value):
        if name in self.__dict__:
            self.__dict__[name] = value
            return
        r, t, v = self.__getparam__(name)
        Py2Xml.val2xml(v, t, value)

    def __setitem__(self, key, value):
        self.__setattr__(key, value)

    def __contains__(self, key):
        try:
            _ = self.__getparam__(key)
        except AttributeError:
            return False
        return True

    def __len__(self):
        return 1

    def __iter__(self):
        self.__iterend__ = False
        return self

    def __next__(self):
        if self.__iterend__:
            raise StopIteration()
        self.__iterend__ = True
        return self

    def update(self, obj):
        if isinstance(obj, Py2Xml):
            d = {k: obj[k] for k in obj.params()}
        else:
            try:
                d = obj.__dict__
            except AttributeError:
                d = obj
            if not isinstance(d, dict):
                raise ValueError('Expected a dictionary or an object with a __dict__ attribute!')

        for k in d:
            try:
                self.__setattr__(k, d[k])
            except AttributeError:
                pass

    def __totree(self):
        ste_root = ET.Element('stexml')
        ste_root.text = '\n'
        ste_root.append(copy.deepcopy(self.__root__))
        ste_root.addprevious(ET.PI('xml-stylesheet', 'type="text/xsl" href="stexml.xsl"'))
        tree = ET.ElementTree(ste_root)
        return tree

    def write(self, filename):
        self.__totree().write(filename, pretty_print=True, encoding='UTF-8', xml_declaration=True)

    def tostring(self):
        return ET.tostring(self.__totree().getroot(), encoding='UTF-8')

    @staticmethod
    def fromstring(string):
        return Xml2Py(ET.fromstring(string).find('object'))

# for backwards-compatiblity (the old Xml2Py class is now obsolete!)
Xml2Py = Py2Xml
