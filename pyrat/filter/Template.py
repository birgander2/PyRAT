import pyrat


# This file shows some templates for the programming of pyrat modules. It contains documented templates for:
#
# - FilterWorker:  Class to process one data set into another one in one single step, using block processing
#                  and multiprocessing. Its functionality is contained in the method 'filter', which is called
#                  on each block of the data.
#
# - WizardWorker:  Class to group / call / run multiple existing pyrat modules to run them in a predefined way.
#                  Used to define entire processing chains in pyrat.
#
# - ImportWorker:  Class for data import
#
# - ExportWorker:  Class for data export
#
# - Worker:        Base class for more general cases. Not so easy to programm, but more flexible. No example
#                  here for the moment, please contact somebody who can explain how to use this one.







# --------------------
# --- FilterWorker ---
# --------------------


class Template(pyrat.FilterWorker):
    """
    Template filterworker class...

    This example does only nonsense, but serves as an example how a PyRat module is programmed.
    These lines are, for example, the documentation of the class, which should be as
    excessive as possible.
    """

    # UNCOMMENT IN NEXT LINE TO MAKE IT APPEAR IN VIEWER!!!
    # gui = {'menu': 'Tools', 'entry': 'Test module'}

    # Defines parameters and their defaults (as well as the dialogbox when calling from the viewer). They
    # will be available later as instance variables (self.win, self.method).
    para = [
        {'var': 'win', 'value': 500, 'type': 'int',  'range': [3, 999], 'text': 'Window size',  'subtext': ['range', 'azimuth']},
        {'var': 'method', 'value': 'original', 'type': 'list', 'range': ['original', 'old', 'new'], 'text': 'Method'}
    ]

    def __init__(self, *args, **kwargs):
        pyrat.FilterWorker.__init__(self, *args, **kwargs)
        self.name = "TEMPLATE"                             # optional, recommended: Name of the class in cli output
        self.blockoverlap = self.win // 2 + 1              # optional: amount of needed overlap between data blocks

        # self.blockprocess = False                        # uncomment to switch off blockprocessing entirely
                                                           # (otherwise data will be always provided in blocks)

        # self.nthreads = 1                                # uncomment to disable multithreading in this class
                                                           # (this is recommended for debugging!)

    # The actual filter: Receives the data in 'array' and returns a filtered version. Metadata can be
    # accessed if necessary (see below). If self.blockprocess is set to True in the constructor, the
    # filter calls are automatically parallelised through blocks in azimuth. If the entire image is needed,
    # don't set self.blockprocess.
    def filter(self, array, *args, **kwargs):                     # the data are provided in 'array'
        meta = kwargs["meta"]                                     # access the meta data dict
        array = self.therealfilter(array, self.win)
        return array

    # Optional definition of a staticmethod: This allows importing the code below for usage outside of PyRat. However,
    # doing everything in filter() is of course also fine. Warning; No access to self here (instance variables)!
    @staticmethod
    def therealfilter(array, win):
        array[:, array.shape[1]/2-win:array.shape[1]/2+win] = 0.0
        return array


# Optional, but highly recommended: definition of a convenience function (avoids appending the run() method)
@pyrat.docstringfrom(Template)                                  # 'steal' docstring from class Template
def template(*args, **kwargs):                                  # use same name as class, but with small letters only
    return Template(*args, **kwargs).run(*args, **kwargs)













# --------------
# --- GROUP ---
# --------------


class MyGroup(pyrat.GroupWorker):
    """
    Template GroupWorker Class. This is basically there to define processing groups
    as pyrat modules, just as to write an independent script.
    """

    # For the usage of the 'gui' and 'para' class variables, see FilterWorker example.

    # The actual method to run: Define a processing strategy of your own choice using python and
    # pyrat modules. In the example, just two pyrat functions are called, one after the other.
    def group(self, *args, **kwargs):
        pyrat.filter.complex2abs()
        pyrat.filter.boxcar(win=self.win)


# Optional, but highly recommended: definition of a convenience function (avoids appending the run() method)
@pyrat.docstringfrom(MyGroup)                                  # 'steal' docstring from class Template
def mygroup(*args, **kwargs):                                  # use same name as class, but with small letters only
    return MyGroup(*args, **kwargs).run(*args, **kwargs)
