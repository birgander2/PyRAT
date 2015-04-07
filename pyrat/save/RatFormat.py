import pyrat
from pyrat.load.tools import RatFile, srat
import logging
import numpy as np


class RatFormat(pyrat.ExportWorker):

    gui = {'menu': 'File', 'entry': 'Save RAT file', 'before': 'Save pixmap'}
    para = [{'var': 'filename', 'value': '', 'type': 'savefile', 'text': 'Save to :'}]

    def __init__(self, *args, **kwargs):
        super(RatFormat, self).__init__(*args, **kwargs)
        self.name = "RAT EXPORT"

    def open(self, *args, **kwargs):
        self.rat = RatFile(self.filename)
        data = pyrat.data.queryLayer(self.layer)
        if isinstance(data, list):
            logging.error("Cannot export multiple layers at once!")
            return False
        else:
            dtype = data['dtype']
            var = self.rat.dtype2var(dtype)
            self.rat.Header.Rat.ndim = data['ndim']
            self.rat.Header.Rat.dim = data['shape'][::-1]
            self.rat.Header.Rat.var = var
            #todo: missing nchannels
            #todo: missing info string
            self.lun = self.rat.write_header()
            return True

    def close(self, *args, **kwargs):
        self.lun.close()
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
            self.lun.seek(1000 + ch * seek_ch + out.itemsize * (kwargs['block'][0]+kwargs['valid'][0]) * nx)
            out[ch, ...].tofile(self.lun)
        return True


def rat(*args, **kwargs):
    return RatFormat(*args, **kwargs).run(**kwargs)
