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


@pyrat.docstringfrom(Coherence)
def coherence(*args, **kwargs):
    return Coherence(*args, **kwargs).run(**kwargs)


class CoherenceMat(pyrat.FilterWorker):
    """
    Calc coherency matrix
    """

    def __init__(self, *args, **kwargs):
        super(CoherenceMat, self).__init__(*args, **kwargs)
        self.name = "CALC COHERENCE MATRIX"
        if 'win' not in self.__dict__:
            self.win = [7, 7]
        self.blockoverlap = self.win[0] / 2 + 1
        self.blockprocess = False
        self.nthreads = 1

    def filter(self, array, *args, **kwargs):

        ntrack = len(array)
        naz, nrg = array[0].shape
        coh = np.empty((ntrack, ntrack, naz, nrg))
        for trk1 in range(ntrack):
            for trk2 in range(ntrack):
                coh[trk1, trk2, ...] = (np.abs(smooth(array[trk1] * np.conj(array[trk2]), self.win) /
                                               np.sqrt(smooth(array[trk1] * np.conj(array[trk1]), self.win) *
                                                       smooth(array[trk2] * np.conj(array[trk2]), self.win))))
        return np.clip(np.nan_to_num(coh), 0.0, 1.0)  # get rid of numerical inaccuracies!


@pyrat.docstringfrom(CoherenceMat)
def coherencemat(*args, **kwargs):
    return CoherenceMat(*args, **kwargs).run(**kwargs)

class CoherenceCov(pyrat.FilterWorker):
    """
    Calc coherency from an inteferometric covariance matrix
    """
    def __init__(self, *args, **kwargs):
        super(CoherenceCov, self).__init__(*args, **kwargs)
        self.name = "CALC InSAR COHERENCE from CM"
        if 'win' not in self.__dict__:
            self.win = [3, 3]
        if 'ch' not in self.__dict__:
            self.ch = [0, 1]
        self.blockoverlap = self.win[0] / 2 + 1
        self.blockprocess = False
        self.nthreads = 1

    def filter(self, array, *args, **kwargs):
        coh = np.abs(smooth(array[self.ch[0],self.ch[1],...], self.win)
                     / np.sqrt(smooth(array[self.ch[0],self.ch[0],...], self.win)
                               * smooth(array[self.ch[1],self.ch[1],...], self.win)))
        return np.clip(np.nan_to_num(coh), 0.0, 1.0)  # get rid of numerical inaccuracies!


@pyrat.docstringfrom(CoherenceCov)
def coherencecov(*args, **kwargs):
    return CoherenceCov(*args, **kwargs).run(**kwargs)

class InterfMat(pyrat.FilterWorker):
    """
    Calc matrix with all interferometers
    """

    def __init__(self, *args, **kwargs):
        super(InterfMat, self).__init__(*args, **kwargs)
        self.name = "CALC INTERFEROMETRY MATRIX"
        if 'win' not in self.__dict__:
            self.win = [7, 7]
        self.blockoverlap = self.win[0] / 2 + 1
        self.blockprocess = False
        self.nthreads = 1

    def filter(self, array, *args, **kwargs):

        ntrack = len(array)
        naz, nrg = array[0].shape
        coh = np.empty((ntrack, ntrack, naz, nrg), dtype=np.complex)
        for trk1 in range(ntrack):
            for trk2 in range(ntrack):
                coh[trk1, trk2, ...] = (smooth(array[trk1] * np.conj(array[trk2]), self.win) /
                                        np.sqrt(smooth(array[trk1] * np.conj(array[trk1]), self.win) *
                                                  smooth(array[trk2] * np.conj(array[trk2]), self.win)))
        return np.nan_to_num(coh)  # get rid of numerical inaccuracies!


@pyrat.docstringfrom(InterfMat)
def interfmat(*args, **kwargs):
    return InterfMat(*args, **kwargs).run(**kwargs)



