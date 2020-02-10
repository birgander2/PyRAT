import shutil

import pyrat
# from pyrat.load.tools import RatFile, srat
from pyrat.lib.ste import RatFile, RatHeader, dtype_dict
import logging
import numpy as np

from pkg_resources import resource_string
from mako.template import Template


class Rat(pyrat.ExportWorker):
    """
    RAT format writer.
    This one should write RAT V2 formats.

    :author: Andreas Reigber
    """
    gui = {'menu': 'File|Export to raster', 'entry': 'RAT (v2)', 'before': 'HDF5 file'}
    para = [{'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save to :', 'extensions': 'RAT (*.rat)'},
            {'var': 'geo_envi_hdr', 'value': False, 'type': 'bool', 'text': 'Write ENVI header'}]

    def __init__(self, *args, **kwargs):
        super(Rat, self).__init__(*args, **kwargs)
        self.name = "RAT EXPORT"
        if len(args) == 1:
            self.file = args[0]

        if 'header' not in self.__dict__:
            self.header = None

    def open(self, *args, **kwargs):
        if isinstance(self.file, tuple):                                       # remove file type if present
            self.file = self.file[0]
        data = pyrat.data.queryLayer(self.layer)
        if data["dtype"] == "bool":                                            # save boolean as uint 8
            data["dtype"] = "uint8"

        if isinstance(data, list):
            logging.error("Cannot export multiple layers at once!")
            return False
        else:
            if self.header is None:
                self.header = RatHeader()
            dtype = data['dtype']
            var = [k for k, v in dtype_dict.items() if v == dtype][0]
            self.header.Rat.ndim = data['ndim']
            self.header.Rat.idl_shape[:data['ndim']] = data['shape'][::-1]
            self.header.Rat.var = var
            # todo: missing nchannels
            # todo: missing info string

            annotation = pyrat.data.getAnnotation(layer=self.layer)
            if "geo_projection" in annotation:
                self.header.Geo.projection = annotation['geo_projection']
            if "geo_min_east" in annotation:
                self.header.Geo.min_east = annotation['geo_min_east']
            if "geo_min_north" in annotation:
                self.header.Geo.min_north = annotation['geo_min_north']
            if "geo_ps_east" in annotation:
                self.header.Geo.ps_east = annotation['geo_ps_east']
            if "geo_ps_north" in annotation:
                self.header.Geo.ps_north = annotation['geo_ps_north']
            if "geo_zone" in annotation:
                self.header.Geo.zone = annotation['geo_zone']

            self.lun = open(self.file, 'wb')
            self.lun.write(self.header)
            self.lun.flush()
            self.rat = RatFile(self.file)
            return True

    def close(self, *args, **kwargs):
        self.lun.close()

        if self.geo_envi_hdr is True:
            logging.info(self.name + '  Writing GEO ENVI Header (.hdr)...')
            hdr = RatFile(self.file).Header
            tmpl = Template(resource_string('pyrat.lib.templates', 'envi_geo_hdr.tpl'))
            envi_hdr = tmpl.render(file=self.file, hdr=hdr)
            with open(self.file+'.hdr','w') as f:
                f.write(envi_hdr)

        return True

    def block_writer(self, array, *args, **kwargs):
        out = array[..., kwargs['valid'][0]:kwargs['valid'][1], :]
        if out.dtype == "bool":                                            # save boolean as uint 8
            out = out.astype("uint8")
        nx = out.shape[-1]
        ny = out.shape[-2]
        nchannels = 1
        for ch in array.shape[:-2]:
            nchannels *= ch
        out = out.reshape((nchannels, ny, nx))
        seek_ch = nx * self.rat.Header.Rat.idl_shape[1] * out.itemsize
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
    gui = {'menu': 'File', 'entry': 'Save PyRAT file', 'before': 'Export to raster'}
    para = [{'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save PyRAT', 'extensions': 'PyRAT (*.ra2)'}]

    def __init__(self, *args, **kwargs):
        super(RatHDF, self).__init__(*args, **kwargs)
        self.name = "PYRAT EXPORT"
        self.nthreads = 1
        if len(args) == 1:
            self.file = args[0]

    def run(self, *args, **kwargs):
        if isinstance(self.file, tuple):                                       # remove file type if present
            self.file = self.file[0]

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
