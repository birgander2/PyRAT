import pyrat
from scipy import misc


class Pixmap(pyrat.ImportWorker):
    # gui = {'menu': 'File', 'entry': 'Import pixmap', 'before': 'File|line1'}
    para = [{'var': 'filename', 'value': '', 'type': 'openfile', 'text': 'file :'}]

    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)
        self.name = 'Pixmap Import'

    def reader(self, *args, **kwargs):
        return misc.imread(self.filename), None


def pixmap(*args, **kwargs):
    Pixmap(*args, **kwargs).run(**kwargs)


class JPG(Pixmap):
    gui = {'menu': 'File|Open pixmap', 'entry': 'JPEG'}


def jpg(*args, **kwargs):
    JPG(*args, **kwargs).run(**kwargs)


class PNG(Pixmap):
    gui = {'menu': 'File|Open pixmap', 'entry': 'PNG'}


def png(*args, **kwargs):
    PNG(*args, **kwargs).run(**kwargs)


class TIFF(Pixmap):
    gui = {'menu': 'File|Open pixmap', 'entry': 'TIFF'}
    key = "TIFF"


def tiff(*args, **kwargs):
    TIFF(*args, **kwargs).run(**kwargs)


