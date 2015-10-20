#!/usr/bin/env python

from __future__ import print_function
from PyQt4 import QtGui, QtCore
QtCore.pyqtRemoveInputHook()
import sys
import logging
import pyrat


def usage():
    logging.basicConfig(format='  %(message)s', level=logging.INFO)
    logging.info(" ")
    logging.info("  Usage: PyRat [OPTION] [RATFILE].")
    logging.info("  Starts the PyRat environment, optinally loading the give filename")
    logging.info(" ")
    logging.info("  --help/-h     : Show this help")
    logging.info("  --batch/-b    : Run in batch mode, not starting the gui")
    logging.info(" ")


def main():
    
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-h', '--help']:
            usage()
            sys.exit()
        elif sys.argv[1] in ['-b', '--batch']:
            # import code
            from pyrat.tools import HistoryConsole
            var = globals().copy()
            var.update(locals())
            shell = HistoryConsole(vars)
            # shell = code.InteractiveConsole(var)
            shell.push("from pyrat import *")
            shell.push("pyrat_init()")
            if len(sys.argv) == 3:
                shell.push("load.rat(filename='"+sys.argv[2]+"')")
            shell.interact("  \n  You are in PyRAT batch mode\n")
            sys.exit()
        else:
            pyrat.pyrat_init()
            pyrat.load.rat(filename=sys.argv[1])
            pyrat.show()
    else:
        pyrat.pyrat_init()
        app = QtGui.QApplication(sys.argv)
        pyrat.app = pyrat.viewer.MainWindow()
        sys.exit(app.exec_())

if __name__ == "__main__":
    main()
