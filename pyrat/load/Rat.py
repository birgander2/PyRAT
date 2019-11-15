import shutil
import tempfile

import pyrat
import logging
import numpy as np
from pyrat.lib.ste import rrat, RatFile


class Rat(pyrat.ImportWorker):
    """
    RAT format reader. Imports the binary data of files in *.rat format into a new PyRAT layer.
    This one should read both RAT V1 and V2 formats.

    :author: Andreas Reigber
    """

    gui = {'menu': 'File|Import raster', 'entry': 'RAT (v1/v2)'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'Open RAT file', 'extensions': 'RAT V1/V2 (*.rat)'}]

    def __init__(self, *args, **kwargs):
        super(Rat, self).__init__(*args, **kwargs)
        self.name = "RAT IMPORT"
        if len(args) == 1:
            self.file = args[0]

    def getmeta(self, *args, **kwargs):
        try:
            file = RatFile(self.file)
            if(file.exists and file.Header.Geo.ps_east != 0.0 and file.Header.Geo.ps_north != 0.0):
                geo = file.Header.Geo
                return {
                        'geo_projection': geo.projection,
                        'geo_min_east': geo.min_east,
                        'geo_min_north': geo.min_north,
                        'geo_ps_east': geo.ps_east,
                        'geo_ps_north': geo.ps_north,
                        'geo_zone': geo.zone,
                        'path': self.file
                    }
        except IOError:
            pass
        return {'path': self.file}

    def getsize(self, *args, **kwargs):
        try:
            file = RatFile(self.file)
            if file.exists:
                size = tuple(file.shape)
                if len(size) == 3 and size[2] < size[0] and size[2] < size[1]:
                    size = (size[2], size[0], size[1])
                if len(size) == 4 and size[2] < size[0] and size[2] < size[1] and size[3] < size[0] \
                        and size[3] < size[1]:
                    size = (size[3], size[2], size[0], size[1])
                del file
                return size
            else:
                logging.error(
                    pyrat.tools.bcolors.FAIL + "IOError: file not found '" + self.file + "'" + pyrat.tools.bcolors.ENDC)
                return False, False
        except IOError:
            logging.error(pyrat.tools.bcolors.FAIL + "IOError: file format not recognised!" + pyrat.tools.bcolors.ENDC)
            return False, False

    def block_reader(self, *args, **kwargs):
        block = kwargs['block']
        array = rrat(self.file, block=block)
        # pixel interleaved vector -> band interleaved vector
        if array.ndim == 3 and array.shape[2] < array.shape[0] and array.shape[2] < array.shape[1]:
            array = np.rollaxis(array, 2)

        # pixel interleaved matrix -> band interleaved vector
        if array.ndim == 4 and array.shape[2] < array.shape[0] and array.shape[2] < array.shape[1] \
                and array.shape[3] < array.shape[0] and array.shape[3] < array.shape[1]:
            array = np.rollaxis(np.rollaxis(array, 2), 3, start=1)
            array = np.reshape(array, (array.shape[0] * array.shape[1], array.shape[2], array.shape[3]))

        array[~np.isfinite(array)] = 0.0

        return array


@pyrat.docstringfrom(Rat)
def rat(*args, **kwargs):
    return Rat(*args, **kwargs).run(*args, **kwargs)


class RatHDF(pyrat.Worker):
    """
    PYRAT format writer  (experimental)

    This is actually a HDF5 file, containing everything from the internal PyRAT structure.
    In particular, meta data is saved in the file, also the preview pyramid of the gui.

    :author: Andreas Reigber
    """
    gui = {'menu': 'File', 'entry': 'Open PyRAT file', 'before': 'Import raster'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'Read PyRAT', 'extensions': 'PyRAT (*.ra2)'}]

    def __init__(self, *args, **kwargs):
        super(RatHDF, self).__init__(*args, **kwargs)
        self.name = "PYRAT IMPORT"
        self.nthreads = 1
        if len(args) == 1:
            self.file = args[0]

    def run(self, *args, **kwargs):
        para = [foo['var'] for foo in self.para]
        self.checkpara(kwargs, para)
        logging.info(
            self.name + '  ' + str(dict((k, v) for k, v in self.__dict__.items() if k in para or k in kwargs)))

        tmp_filename = tempfile.mktemp(suffix='.hd5', prefix='pyrat_', dir=pyrat.data.tmpdir)
        shutil.copyfile(self.file, tmp_filename)
        lay = pyrat.data.addLayer(file=tmp_filename)
        pyrat.activate(lay)
        return lay

@pyrat.docstringfrom(RatHDF)
def rathdf(*args, **kwargs):
    return RatHDF(*args, **kwargs).run(*args, **kwargs)
