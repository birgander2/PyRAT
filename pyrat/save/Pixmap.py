import pyrat
from pyrat.viewer.tools import sarscale, phascale, cohscale
import logging
import numpy as np
from scipy import misc
from pyrat.tools import colortables


class Pixmap(pyrat.ExportWorker):
    """
    Export to various pixmap formats.

    For possible parameters see source code ;-)
    """
    para = [
        {'var': 'file', 'value': '', 'type': 'savefile', 'text': 'Save to :'},
        {'var': 'chscl', 'value': True, 'type': 'bool', 'text': 'Scale channels indiviually'},
        {'var': 'method', 'value': 'amplitude', 'type': 'list',
         'range': ['amplitude', 'intensity', 'phase', 'coherence'], 'text': 'Method'},
        {'var': 'scaling', 'value': 2.5, 'type': 'float', 'range': [0.1, 20.0], 'text': 'SAR scaling factor'},
        {'var': 'palette', 'value': 'bw linear', 'type': 'list', 'range': colortables()[0], 'text': 'Color table'}
    ]
    key = None

    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)
        self.name = "EXPORT TO PIXMAP"
        if len(args) == 1:
            self.file = args[0]
        if 'order' not in self.__dict__:
            self.order = [0, 1, 2]

    def writer(self, array, *args, **kwargs):
        if isinstance(self.file, tuple):                                       # remove file type if present
            self.file = self.file[0]

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
            out = colortables(self.palette)[1][out]

        try:
            misc.imsave(self.file, out, format=self.key)
            return True
        except IOError as err:
            logging.error("ERROR:" + str(err))
            return False
        else:
            logging.error("UNKNOWN ERROR")
            return False


@pyrat.docstringfrom(Pixmap)
def pixmap(*args, **kwargs):
    Pixmap(*args, **kwargs).run(*args, **kwargs)


class JPG(Pixmap):
    """
    JPG format writer
    """
    gui = {'menu': 'File|Save as pixmap', 'entry': 'JPEG'}
    key = "JPEG"


@pyrat.docstringfrom(JPG)
def jpg(*args, **kwargs):
    JPG(*args, **kwargs).run(*args, **kwargs)


class PNG(Pixmap):
    """
    PNG format writer
    """
    gui = {'menu': 'File|Save as pixmap', 'entry': 'PNG'}
    key = "PNG"


@pyrat.docstringfrom(PNG)
def png(*args, **kwargs):
    PNG(*args, **kwargs).run(*args, **kwargs)


class TIFF(Pixmap):
    """
    TIFF format writer
    """
    gui = {'menu': 'File|Save as pixmap', 'entry': 'TIFF'}
    key = "TIFF"


@pyrat.docstringfrom(TIFF)
def tiff(*args, **kwargs):
    TIFF(*args, **kwargs).run(*args, **kwargs)


class PDF(Pixmap):
    """
    PDF format writer
    """
    gui = {'menu': 'File|Save as pixmap', 'entry': 'PDF'}
    key = "PDF"


@pyrat.docstringfrom(PDF)
def pdf(*args, **kwargs):
    PDF(*args, **kwargs).run(*args, **kwargs)


class EPS(Pixmap):
    """
    EPS format writer
    """
    gui = {'menu': 'File|Save as pixmap', 'entry': 'EPS'}
    key = "EPS"


@pyrat.docstringfrom(EPS)
def eps(*args, **kwargs):
    EPS(*args, **kwargs).run(*args, **kwargs)
