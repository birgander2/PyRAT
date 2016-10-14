import shutil

import pyrat
from pyrat.load.tools import RatFile, srat
import logging
import numpy as np

from  pkg_resources import resource_string
from mako.template import Template


class Rat(pyrat.ExportWorker):
    """
    RAT format writer.
    This one should write RAT V2 formats.

    :author: Andreas Reigber
    """
    gui = {'menu': 'File', 'entry': 'Save RAT file', 'before': 'Save pixmap'}
    para = [{'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save to :'}]

    def __init__(self, *args, **kwargs):
        super(Rat, self).__init__(*args, **kwargs)
        self.name = "RAT EXPORT"
        if len(args) == 1:
            self.file = args[0]

        if 'header' not in self.__dict__:
            self.header = None

    def open(self, *args, **kwargs):
        self.rat = RatFile(self.file)
        data = pyrat.data.queryLayer(self.layer)
        if isinstance(data, list):
            logging.error("Cannot export multiple layers at once!")
            return False
        else:
            dtype = data['dtype']
            var = self.rat.dtype2var(dtype)
            if self.header is not None:
                self.rat.Header = self.header
            self.rat.Header.Rat.ndim = data['ndim']
            self.rat.Header.Rat.dim[:data['ndim']] = data['shape'][::-1]
            self.rat.Header.Rat.var = var
            # todo: missing nchannels
            # todo: missing info string
            self.lun = self.rat.write_header()
            return True

    def close(self, *args, **kwargs):
        self.lun.close()

        if 'geo_envi_hdr' in self.__dict__ and self.geo_envi_hdr:
            logging.info(self.name + '  Writing GEO ENVI Header (.hdr)...')
            hdr = RatFile(self.file).Header
            tmpl = Template(resource_string('pyrat.templates', 'envi_geo_hdr.tpl'))
            envi_hdr = tmpl.render(file=self.file, hdr=hdr)
            with open(self.file+'.hdr','w') as f:
                f.write(envi_hdr)

        return True

    def block_writer(self, array, *args, **kwargs):
        out = array[..., kwargs['valid'][0]:kwargs['valid'][1], :]

        nx = out.shape[-1]
        ny = out.shape[-2]
        nchannels = 1
        for ch in array.shape[:-2]:
            nchannels *= ch
        out = out.reshape((nchannels, ny, nx))
        seek_ch = nx * self.rat.Header.Rat.dim[1] * out.itemsize
        for ch in range(nchannels):
            self.lun.seek(1000 + ch * seek_ch + out.itemsize * (kwargs['block'][0] + kwargs['valid'][0]) * nx)
            out[ch, ...].tofile(self.lun)
        return True


@pyrat.docstringfrom(Rat)
def rat(*args, **kwargs):
    return Rat(*args, **kwargs).run(*args, **kwargs)


class RatHDF(pyrat.Worker):
    """
    PYRAT format writer (experimental)

    This is actually a HDF5 file, containing everything from the internal PyRAT structure.
    In particular, meta data is saved in the file, also the preview pyramid of the gui.

    :author: Andreas Reigber
    """
    gui = {'menu': 'File', 'entry': 'Save PYRAT file', 'before': 'Save pixmap'}
    para = [{'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save PyRAT'}]

    def __init__(self, *args, **kwargs):
        super(RatHDF, self).__init__(*args, **kwargs)
        self.name = "PYRAT EXPORT"
        self.nthreads = 1
        if len(args) == 1:
            self.file = args[0]

    def run(self, *args, **kwargs):
        para = [foo['var'] for foo in self.para]
        self.checkpara(kwargs, para)
        logging.info(
            self.name + '  ' + str(dict((k, v) for k, v in self.__dict__.items() if k in para or k in kwargs)))

        if isinstance(self.layer, list):
            logging.error("Cannot export (yet) multiple layers at once!")
            return
        if pyrat.data.layers[self.layer].attrs['_type'] == 'Memory':
            logging.error("ERROR: Cannot save (yet) memory layers!")
            return
        pyrat.data.layers[self.layer].file.flush()
        shutil.copyfile(pyrat.data.layers[self.layer].fn, self.file)


@pyrat.docstringfrom(RatHDF)
def rathdf(*args, **kwargs):
    return RatHDF(*args, **kwargs).run(*args, **kwargs)
