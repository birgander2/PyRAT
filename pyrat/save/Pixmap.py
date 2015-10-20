import pyrat
from pyrat.viewer.tools import sarscale, phascale, cohscale
import logging
import numpy as np
from scipy import misc
from pyrat.tools import colortables


class Pixmap(pyrat.ExportWorker):
    para = [
        {'var': 'filename', 'value': '', 'type': 'savefile', 'text': 'Save to :'},
        {'var': 'chscl', 'value': True, 'type': 'bool', 'text': 'Scale channels indiviually'},
        {'var': 'method', 'value': 'amplitude', 'type': 'list', 'range': ['amplitude', 'intensity', 'phase', 'coherence'], 'text': 'Method'},
        {'var': 'scaling', 'value': 2.5, 'type': 'float', 'range': [0.1, 20.0], 'text': 'SAR scaling factor'},
        {'var': 'palette', 'value': 'bw linear', 'type': 'list', 'range': colortables()[0], 'text': 'Color table'}
    ]
    key = None

    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)
        self.name = "EXPORT TO PIXMAP"
        if 'order' not in self.__dict__:
            self.order = [0, 1, 2]

    def writer(self, array, *args, **kwargs):
        if self.method == 'amplitude' or self.method == 'intensity' and np.iscomplexobj(array):
            array = np.abs(array)

        array = np.squeeze(array)
        out = np.zeros_like(array, dtype='uint8')
        if array.ndim == 3 and self.chscl == True:
            for k in range(array.shape[0]):
                array[k, ...] = array[k, ...] / np.mean(array[k, ...])

        if self.method == 'intensity' or self.method == 'amplitude':
            if self.method == 'amplitude':
                array **= 0.7
            if self.method == 'intensity':
                array **= 0.35
            print(self.scaling)
            out = sarscale(array, factor=self.scaling)
        elif self.method == 'phase':
            out = phascale(array)
        elif self.method == 'coherence':
            out = cohscale(array)
        else:
            logging.error("Scaling method unknown")
            return False

        if array.ndim == 3:
            out = out[self.order, ...]
        else:
            if self.palette != 'bw linear':                                   # apply color palette
                out = colortables(self.palette)[1][out]

        try:
            misc.imsave(self.filename, out, format=self.key)
            return True
        except IOError as err:
            logging.error("ERROR:"+str(err))
            return False
        else:
            logging.error("UNKNOWN ERROR")
            return False

def pixmap(*args, **kwargs):
    Pixmap(*args, **kwargs).run(**kwargs)


class JPG(Pixmap):
    gui = {'menu': 'File|Save pixmap', 'entry': 'JPEG'}
    key = "JPEG"


def jpg(*args, **kwargs):
    JPG(*args, **kwargs).run(**kwargs)


class PNG(Pixmap):
    gui = {'menu': 'File|Save pixmap', 'entry': 'PNG'}
    key = "PNG"


def png(*args, **kwargs):
    PNG(*args, **kwargs).run(**kwargs)


class TIFF(Pixmap):
    gui = {'menu': 'File|Save pixmap', 'entry': 'TIFF'}
    key = "TIFF"


def tiff(*args, **kwargs):
    TIFF(*args, **kwargs).run(**kwargs)


class PDF(Pixmap):
    gui = {'menu': 'File|Save pixmap', 'entry': 'PDF'}
    key = "PDF"


def pdf(*args, **kwargs):
    PDF(*args, **kwargs).run(**kwargs)


class EPS(Pixmap):
    gui = {'menu': 'File|Save pixmap', 'entry': 'EPS'}
    key = "EPS"


def eps(*args, **kwargs):
    EPS(*args, **kwargs).run(**kwargs)
