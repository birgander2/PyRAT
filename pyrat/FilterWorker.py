from __future__ import print_function
import pyrat
import logging
import copy
from pyrat.tools import ProgressBar


def exec_parallel(args):
    return args[0].filter(args[1], meta=args[2], block=args[3]), args[2]


class FilterWorker(pyrat.Worker):
    def __init__(self, *args, **kwargs):
        super(FilterWorker, self).__init__(*args, **kwargs)

    def run(self, *args, **kwargs):
        """
        Main routine doing the (parallel) block processing, calling the (overloaded) filter method.
        Before the actual processing, pre() is called, as well as post() after completing it.
        """
        para = [foo['var'] for foo in self.para]
        self.checkpara(kwargs, para)
        logging.info(self.name + '  ' + str(dict((k, v) for k, v in self.__dict__.items()
                                                 if k in para or k in kwargs)))
        if self.checkinput():
            self.pre(*args, **kwargs)
            newlayer = self.layer_process(self.filter, silent=False, **kwargs)
            pyrat.data.activateLayer(newlayer)
            self.post(*args, **kwargs)
            return newlayer
        else:
            return False

    def filter(self, data, *args, **kwargs):
        """
        The actual filter routine (to be overloaded)
        """
        logging.error(self.name + ': No filter method defined')
        return False

    def pre(self, *args, **kwargs):
        """
        The preprocessing routine (to be overloaded)
        """
        pass

    def post(self, *args, **kwargs):
        """
        The postprocessing routine (to be overloaded)
        """
        pass

    def info(self):
        print(self.__doc__)
