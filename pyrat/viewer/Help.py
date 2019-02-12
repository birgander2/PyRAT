from PyQt5 import QtWidgets
import pyrat
import sys


class About(pyrat.Worker):
    gui = {'menu': 'Help', 'entry': 'About PyRat'}

    @classmethod
    def guirun(cls, viewer):
        message = """
        PyRAT - Radar Tools Collection
        ---
        (c) 2015-2019 by the PyRAT development team
        ---
        PyRAT version: v%s
        running on: %s / python %s
        ---
        """ % (str(pyrat.__version__), sys.platform, sys.version[0:3])
        foo = QtWidgets.QMessageBox(parent=viewer)
        foo.setIcon(1)
        foo.setText(message)
        foo.exec_()

