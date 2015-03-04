import pyrat
from scipy import misc


class Pixmap(pyrat.ImportWorker):
    gui = {'menu': 'File', 'entry': 'Open pixmap', 'before': 'File|line1'}
    para = [{'var': 'filename', 'value': '', 'type': 'openfile', 'text': 'file :'}]

    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)
        self.name = 'Pixmap Import'

    def reader(self, *args, **kwargs):
        return misc.imread(self.filename), None


def pixmap(*args, **kwargs):
    Pixmap(*args, **kwargs).run(**kwargs)
