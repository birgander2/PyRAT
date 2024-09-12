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
        (c) 2015-2024 by the PyRAT development team
        ---
        PyRAT version: v%s
        running on: %s / python %s
        ---
        """ % (str(pyrat.__version__), sys.platform, str(sys.version_info.major) + "." + str(sys.version_info.minor))
        foo = QtWidgets.QMessageBox(parent=viewer)
        foo.setIcon(1)
        foo.setText(message)
        foo.exec_()

