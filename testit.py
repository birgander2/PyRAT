#!/usr/bin/env python

import PyRat

import os
os.chdir('/Users/mneumann/sar/py/pyrat')


PyRat.init()
l1 = PyRat.Import.Rat(filename="alling_k3l.rat").run()
PyRat.Data.setAnnotation({'CH_pol': ['HH','VV','XX']})
l2 = PyRat.Filter.Lex2Pauli().run()
l2 = PyRat.Filter.Vec2Mat().run(memory=True)
l3 = PyRat.Filter.Despeckle.Boxcar().run()
l4 = PyRat.Filter.Eigen().run()
l5 = PyRat.Filter.Entalpani().run()
PyRat.Export.Rat(filename="alling_pyout.rat").run()
PyRat.Data.close()


#PyRat.Viewer.MainWindow()