#!/usr/bin/python
from PyQt4 import QtGui, QtCore
QtCore.pyqtRemoveInputHook()
import sys
import PyRat

import logging
logging.basicConfig(level=logging.DEBUG)

def usage():
    print()
    print("PyRAT - Radar Tools (v0.1)")
    print("==========================")
    print("Command line arguments:")
    print("--help/-h  : Show some help")
    print("--nogui    : Run in batch mode (not working yet")
    print()

def main():

    if len(sys.argv) > 1:
        if sys.argv[1] == '-h' or sys.argv[1] == '--help' :
            usage()
            sys.exit()
        if sys.argv[1] == '--nogui' :
            print("Now batch mode should start up...")
            print("Sorry, not implemented yet")
            sys.exit()

    app   = QtGui.QApplication(sys.argv)
    PyRat.init()
    print("PyRat.init()")
    pyrat = PyRat.Viewer.MainWindow()
    print("pyrat = PyRat.Viewer.MainWindow()")
    sys.exit(app.exec_())
    print("sys.exit(app.exec_())")

if __name__ == "__main__":
    main()
