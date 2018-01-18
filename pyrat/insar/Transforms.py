import pyrat
import numpy as np


class Phase(pyrat.FilterWorker):
    """
    Calc InSAR phase between two images, or extract it from complex interferogram
    """
    gui = {'menu': 'InSAR|Transform', 'entry': 'Extract phase'}

    def __init__(self, *args, **kwargs):
        super(Phase, self).__init__(*args, **kwargs)
        self.name = "Extract InSAR PHASE"

    def filter(self, array, *args, **kwargs):
        if not np.iscomplexobj(array):
            raise ValueError("Data not complex")
        if isinstance(array, list):
            return np.angle(array[0] * np.conj(array[1]))
        else:
            return np.angle(array)


@pyrat.docstringfrom(Phase)
def phase(*args, **kwargs):
    return Phase(*args, **kwargs).run(*args, **kwargs)
