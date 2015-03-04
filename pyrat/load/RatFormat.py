import pyrat
import logging
import numpy as np
from PyQt4 import QtCore, QtGui
from pyrat.load.tools import rrat, RatFile


class RatFormat(pyrat.ImportWorker):

    gui = {'menu': 'File', 'entry': 'Open RAT file', 'before': 'Open external'}
    para = [{'var': 'filename', 'value': '', 'type': 'openfile', 'text': ''}]

    def __init__(self, *args, **kwargs):
        super(RatFormat, self).__init__(*args, **kwargs)
        self.name = "RAT IMPORT"

    def getsize(self, *args, **kwargs):
        try:
            file = RatFile(self.filename)
            size = tuple(file.dim[:file.ndim][::-1])
            if len(size) == 3 and size[2] < size[0] and size[2] < size[1]:
                size = (size[2], size[0], size[1])
            if len(size) == 4 and size[2] < size[0] and size[2] < size[1] and size[3] < size[0] and size[3] < size[1]:
                size = (size[2], size[3], size[0], size[1])
            del file
            return size
        except IOError:
            logging.error("IOError: file format not recognised!")
            return False

    def block_reader(self, *args, **kwargs):
        block = kwargs['block']
        block = [block[2], block[0], block[3]-block[2], block[1]-block[0]]

        array = rrat(self.filename, block=block)

        # pixel interleaved vector -> band interleaved vector
        if array.ndim == 3 and array.shape[2] < array.shape[0] and array.shape[2] < array.shape[1]:
            array = np.rollaxis(array, 2)

        # pixel interleaved matrix -> band interleaved vector
        if array.ndim == 4 and array.shape[2] < array.shape[0] and array.shape[2] < array.shape[1]  \
                and array.shape[3] < array.shape[0] and array.shape[3] < array.shape[1]:
            array = np.rollaxis(np.rollaxis(array, 2), 3, start=1)
            array = np.reshape(array, (array.shape[0] * array.shape[1], array.shape[2], array.shape[3]))
        array[~np.isfinite(array)] = 0.0

        return array


def rat(*args, **kwargs):
    return RatFormat(*args, **kwargs).run(**kwargs)
