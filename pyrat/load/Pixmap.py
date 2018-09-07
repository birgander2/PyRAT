import pyrat
from PIL import Image
import numpy as np


class Pixmap(pyrat.ImportWorker):
    """
    Reader for varoious graphic formats (all what is supported by pillow)
    """
    # gui = {'menu': 'File', 'entry': 'Import pixmap', 'before': 'File|line1'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'file :'}]

    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)
        self.name = 'Pixmap Import'
        self.scaling_hint = 'min->max'
        if len(args) == 1:
            self.file = args[0]

    def reader(self, *args, **kwargs):
        data = np.array(Image.open(self.file))
        if data.ndim == 3:   # color image
            data = np.rollaxis(data, 2, start=0)
        return data, None


@pyrat.docstringfrom(Pixmap)
def pixmap(*args, **kwargs):
    Pixmap(*args, **kwargs).run(*args, **kwargs)


class JPG(Pixmap):
    """
    JPG format reader
    """
    gui = {'menu': 'File|Import pixmap', 'entry': 'JPEG'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'file :', 'extensions': 'JPEG file (*.jpg, *.jpeg)'}]


@pyrat.docstringfrom(JPG)
def jpg(*args, **kwargs):
    JPG(*args, **kwargs).run(*args, **kwargs)


class PNG(Pixmap):
    """
    PNG format reader
    """
    gui = {'menu': 'File|Import pixmap', 'entry': 'PNG'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'file :', 'extensions': 'PNG file (*.png)'}]


@pyrat.docstringfrom(PNG)
def png(*args, **kwargs):
    PNG(*args, **kwargs).run(*args, **kwargs)


class TIFF(Pixmap):
    """
    Tiff format reader
    """
    gui = {'menu': 'File|Import pixmap', 'entry': 'TIFF'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'file :', 'extensions': 'TIFF file (*.tif *.tiff)'}]
    key = "TIFF"


@pyrat.docstringfrom(TIFF)
def tiff(*args, **kwargs):
    TIFF(*args, **kwargs).run(*args, **kwargs)
