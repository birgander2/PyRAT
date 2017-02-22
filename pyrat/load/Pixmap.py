import pyrat
from scipy import misc


class Pixmap(pyrat.ImportWorker):
    """
    Reader for varoious graphic formats (all what is supported by pillow)
    """
    # gui = {'menu': 'File', 'entry': 'Import pixmap', 'before': 'File|line1'}
    para = [{'var': 'file', 'value': '', 'type': 'openfile', 'text': 'file :'}]

    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)
        self.name = 'Pixmap Import'
        if len(args) == 1:
            self.file = args[0]

    def reader(self, *args, **kwargs):
        return misc.imread(self.file), None


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
