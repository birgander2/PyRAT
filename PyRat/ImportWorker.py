import PyRat
import logging, pdb
from yapsy.IPlugin import IPlugin
from PyQt4 import QtGui, QtCore

class ImportWorker(IPlugin):
    def __init__(self, *args, **kwargs):
        for (k, v) in kwargs.items():           # copy keywords to self
            setattr(self, k, v)
        self.name = 'UNKNOWN'
        
    def run(self, *args, **kwargs):
        logging.info(self.name + '  '+str(dict((k, v) for k,v in self.__dict__.items() if k not in ['name','blockprocessing'])))
        track = None
        annotation = {}
        newlayer = PyRat.Data.addLayer(self.reader(track=track, attrs=annotation, *args, **kwargs))
        PyRat.Data.setAnnotation(annotation, layers=newlayer)
        PyRat.Data.activateLayer(newlayer)
        return newlayer
        
    def reader(self, filename, *args, **kwargs):
        return False

    def registerGUI(self, viewer):
        action = QtGui.QAction(viewer.plugin.entry, viewer, shortcut=viewer.plugin.shortcut)
        action.setStatusTip(viewer.plugin.tooltip)
        viewer.connect(action, QtCore.SIGNAL('triggered()'), lambda: self.guirun(viewer))
        viewer.menue[viewer.plugin.menu].insertAction(viewer.exitAct,action)  # insert entry before ...
    
    def guirun(self, viewer):
        self.filename = str(QtGui.QFileDialog(viewer).getOpenFileName())
        viewer.statusBar.setMessage(message="Reading "+self.filename, colour = 'R')
        layer = self.run()
        viewer.genPyramid(PyRat.Data.getData(layer=layer))
        viewer.processPicture(fitwin=1)
        viewer.statusBar.setMessage(colour = 'G')
