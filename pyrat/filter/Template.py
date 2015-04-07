from __future__ import print_function
import pyrat


class Template(pyrat.FilterWorker):
    """
    Template filter class...

    It does only nonsense, but serves as an example how a PyRat module is programmed.
    These lines are, for example, the documentation of the class, which should be as
    excessive as possible.
    """

    # Optional line, adds the class to the viewer menu.

    # COMMENT IN NEXT LINE TO MAKE IT APPEAR IN VIEWER!!!
    # gui = {'menu': 'Tools', 'entry': 'Test module'}

    # Defines parameters and their defaults (as well as the dialogbox when calling from the viewer)
    para = [
        {'var': 'win', 'value': 500, 'type': 'int',  'range': [3, 999], 'text': 'Window size',  'subtext': ['range', 'azimuth']},
        {'var': 'method', 'value': 'original', 'type': 'list', 'range': ['original', 'old', 'new'], 'text': 'Method'}
    ]

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.blockprocess = True

    # The actual filter: Receives the data in 'array' and returns a filtered version. Metadata can be
    # accessed if necessary (see below). If self.blockprocess is set to True in the constructor, the
    # filter calls are automatically parallelised through blocks in azimuth. If the entire image is needed,
    # don't set self.blockprocess.
    def filter(self, array, *args, **kwargs):
        meta = kwargs["meta"]
        array = self.therealfilter(array, self.win)
        return array

    # Optional definition of a staticmethod: This allows importing the code below for usage outside of PyRat. However,
    # doing everything in filter() is of course also fine. Warning; No access to self here (instance variables)!
    @staticmethod
    def therealfilter(array, win):
        array[:, array.shape[1]/2-win:array.shape[1]/2+win] = 0.0
        return array


    # Optional: definition of a convenience function (avoids appending the run() method)
def template(*args, **kwargs):
    return Template(*args, **kwargs).run(**kwargs)
