from __future__ import print_function
import pyrat
import logging
from PyQt4 import QtGui, QtCore


class WizardWorker(pyrat.Worker):
    def __init__(self, *args, **kwargs):
        super(WizardWorker, self).__init__(*args, **kwargs)
        for (k, v) in self.para.items():  # copy keywords to self
            setattr(self, k, v['value'])
        for (k, v) in kwargs.items():  # copy keywords to self
            setattr(self, k, v)
        self.name = self.__class__.__name__

    def run(self, *args, **kwargs):
        logging.info(
            self.name + '  ' + str(dict((k, v) for k, v in self.__dict__.items() if k in self.para or k in kwargs)))
        self.wizard(*args, **kwargs)
        newlayer = pyrat.data.active
        if len(newlayer) == 1:
            newlayer = newlayer[0]
        return newlayer

    def wizard(self, *args, **kwargs):
        print("ERROR: WizardWorker.wizard() not overloaded")
        return False, False
