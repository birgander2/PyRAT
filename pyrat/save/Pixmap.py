import pyrat
from pyrat.viewer.tools import sarscale, phascale, cohscale
import logging
import numpy as np
from scipy import misc


class Pixmap(pyrat.ExportWorker):
    gui = {'menu': 'File', 'entry': 'Save pixmap', 'before': 'Exit'}
    para = [
        {'var': 'filename', 'value': '', 'type': 'savefile', 'text': 'Save to :'},
        {'var': 'chscl', 'value': True, 'type': 'bool', 'text': 'Scale channels indiviually'},
        {'var': 'scale', 'value': 'SAR', 'type': 'list', 'range': ['SAR', 'phase', 'coherence'], 'text': 'Scaling'}
    ]

    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)
        self.name = "EXPORT TO PIXMAP"
        if 'scale' not in self.__dict__:
            self.scale = 'SAR'
        if 'chscl' not in self.__dict__:
            self.chscl = True
        if 'order' not in self.__dict__:
            self.order = [0, 1, 2]

    def writer(self, array, *args, **kwargs):
        if self.scale == 'SAR' and np.iscomplexobj(array):
            array = np.abs(array)

        array = np.squeeze(array)
        out = np.zeros_like(array, dtype='uint8')
        if array.ndim == 3 and self.chscl == True:
            for k in range(array.shape[0]):
                array[k, ...] = array[k, ...] / np.mean(array[k, ...])

        if self.scale == 'SAR':
            out = sarscale(array)
        elif self.scale == 'phase':
            out = phascale(array)
        elif self.scale == 'coherence':
            out = cohscale(array)
        else:
            logging.error("Scaling method unknown")
            return False

        if array.ndim == 3:
            out = out[self.order, ...]

        misc.imsave(self.filename, out)

        return True


def pixmap(*args, **kwargs):
    Pixmap(*args, **kwargs).run(**kwargs)
