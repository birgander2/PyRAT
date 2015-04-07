from __future__ import print_function
import pyrat
from pyrat.filter.tools import rebin


class Presum(pyrat.FilterWorker):
    """
    Presumming of data arrays.

    :author: Andreas Reigber
    """
    gui = {'menu': 'Tools', 'entry': 'Presumming'}
    para = [
        {'var': 'subx', 'value': 4, 'type': 'int',  'range': [0, 999], 'text': 'Presumming range'},
        {'var': 'suby', 'value': 4, 'type': 'int',  'range': [0, 999], 'text': 'Presumming azimuth'}
        ]

    def __init__(self, *args, **kwargs):
        super(Presum, self).__init__(*args, **kwargs)
        self.name = "PRESUMMING"
        self.blockprocess = False

    def filter(self, array, *args, **kwargs):

        odim = list(array.shape)
        odim[-2] //= self.suby
        odim[-1] //= self.subx

        array = rebin(array[..., 0:odim[-2] * self.suby, 0:odim[-1] * self.subx], odim)

        return array

def presum(*args, **kwargs):
    return Presum(*args, **kwargs).run(**kwargs)


# class PresumPlus(pyrat.FilterWorker):
#     """
#     Presumming of data arrays.
#
#     :author: Andreas Reigber
#     """
#     gui = {'menu': 'Tools', 'entry': 'Presumming++'}
#     para = [
#         {'var': 'subx', 'value': 4, 'type': 'int',  'range': [0, 999], 'text': 'Presumming range'},
#         {'var': 'suby', 'value': 4, 'type': 'int',  'range': [0, 999], 'text': 'Presumming azimuth'}
#         ]
#
#     def __init__(self, *args, **kwargs):
#         super(PresumPlus, self).__init__(*args, **kwargs)
#         self.name = "PRESUMMING"
#         self.blockprocess = True
#         self.nthreads = 1
#
#         self.blockmul = self.suby
#         self.outshape = list(pyrat.data.shape)
#         self.outshape[-2] //= self.suby
#         self.outshape[-1] //= self.subx
#
#
#         stop()
#
#     def filter(self, array, *args, **kwargs):
#
#         odim = list(array.shape)
#         odim[-2] //= self.suby
#         odim[-1] //= self.subx
#
#         array = rebin(array[..., 0:odim[-2] * self.suby, 0:odim[-1] * self.subx], odim)
#
#         return array
