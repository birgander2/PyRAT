import pyrat
import numpy as np
from pyrat.filter.tools import smooth


class Coherence(pyrat.FilterWorker):
    """
    Calc InSAR phase between two images
    """

    def __init__(self, *args, **kwargs):
        super(Coherence, self).__init__(*args, **kwargs)
        self.name = "CALC InSAR COHERENCE"
        if 'win' not in self.__dict__:
            self.win = [7, 7]
        self.blockoverlap = self.win[0] / 2 + 1
        self.blockprocess = True
        self.nthreads = 1

    def filter(self, array, *args, **kwargs):
        coh = np.abs(smooth(array[0] * np.conj(array[1]), self.win)
                     / np.sqrt(smooth(array[0] * np.conj(array[0]), self.win)
                               * smooth(array[1] * np.conj(array[1]), self.win)))
        return np.clip(np.nan_to_num(coh), 0.0, 1.0)  # get rid of numerical inaccuracies!
