from __future__ import print_function
from PyQt4 import QtCore, QtGui
import pyrat
import sys


class About(pyrat.Worker):
    gui = {'menu': 'Help', 'entry': 'About PyRat'}

    @classmethod
    def guirun(cls, viewer):
        message = """
        PyRat - Radar Tools Collection
        ---
        (c) 2014-2015 by the RAT development team
        ---
        PyRat version: v%s
        running on: %s / python %s
        ---
        """ % (str(pyrat.__version__), sys.platform, sys.version[0:3])
        foo = QtGui.QMessageBox(parent=viewer)
        foo.setIcon(1)
        foo.setText(message)
        foo.exec_()

