from __future__ import print_function

import pyrat
import numpy as np

try:
    from pyrat.lib.nlsar import sarnlsar
    haveLib = True
except ImportError:
    haveLib = False

class NLPolSAR(pyrat.FilterWorker):
    """
    Non-local means filter for PolSAR data
    (wrapper only: C implementation by Charles Deledalle under lib/nlsar/).

    :author: Marc Jaeger
    :param array: The covariance matrix to filter (4D np.ndarray), optionally also a second activated layer that represents a region of pure speckle noise (4D np.ndarray) used to derive a more accurate statistical signal model.
    :type array: complex64
    :param enl: The effective number of looks of the input (default: 1)
    :type enl: integer
    :param hW: The radius within which to search for matching patches (default 12)
    :type hW: int
    :param hP: The radius of each patch (default 5)
    :type hP: int
    :returns: filtered covariance matrix (4D np.ndarray) and the estimated effective number of looks (2D np.ndarray)
    """

    gui = {'menu': 'PolSAR|Speckle filter', 'entry': 'Non-local Means'}
    para = [
        {'var': 'enl', 'value': 1, 'type': 'int', 'range': [1, 999], 'text': 'Effective # of looks'},
        {'var': 'hW', 'value': 12, 'type': 'int', 'range': [1, 999], 'text': 'Search radius'},
        {'var': 'hP', 'value': 5, 'type': 'int', 'range': [1, 999], 'text': 'Patch radius'}
    ]

    def __init__(self, *args, **kwargs):
        super(NLPolSAR, self).__init__(*args, **kwargs)
        self.name = "NON-LOCAL MEANS"
        self.blockprocess = False

        if 'enl' not in self.__dict__:
            self.enl = 1
        if 'hW' not in self.__dict__:
            self.hW = 12
        if 'hP' not in self.__dict__:
            self.hP = 5


    def filter(self, array, *args, **kwargs):

        if not haveLib:
            raise Exception('NLSAR C library not found!\n' + \
                            'Please compile it by issuing the command\n' + \
                            '"./configure --disable-matlab && make" in in the "lib/nlsar"\n' + \
                            'directory of the PyRAT distribution.')

        try:
            if isinstance(array,list):
                if array[0].ndim == 2:
                    array[0] = array[0][np.newaxis,np.newaxis,...]
                result = sarnlsar(array[0].T, self.enl, 1, self.hW, self.hP, array[1].T)
            else:
                if array.ndim == 2:
                    array = array[np.newaxis,np.newaxis,...]
                result = sarnlsar(array.T, self.enl, 1, self.hW, self.hP)
        except Warning:
            print(str(Warning))

        if result[0].ndim == 2:
            result[0] = result[0][...,np.newaxis,np.newaxis]

        return [result[0].T,result[1].T]
