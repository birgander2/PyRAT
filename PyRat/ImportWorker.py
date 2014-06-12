import PyRat
import logging
import numpy as np
import pdb
from yapsy.IPlugin import IPlugin
from PyQt4 import QtGui, QtCore

class ImportWorker(IPlugin):
    def __init__(self, *args, **kwargs):
        for (k, v) in kwargs.items():           # copy keywords to self
            setattr(self, k, v)
        self.name = 'UNKNOWN'
        self.block = 'D'

    def run(self, *args, **kwargs):
        logging.info(self.name + '  '+str(dict((k, v) for k,v in self.__dict__.items() if k not in ['name','block'])))
        data, meta = self.reader(*args, **kwargs)

        if data == None and meta == None:
            logging.warning("Nothing imported!!!!!!!!!!")
            return None
        else:
            if data == None: logging.warning("No image data imported")
            if meta == None: logging.warning("No meta data imported")
            newlayer = PyRat.Data.addLayer(data, block=self.block, *args, **kwargs)
            PyRat.Data.setAnnotation(meta, layer=newlayer)
            if "name" in self.__dict__:
                PyRat.Data.setAnnotation({"Created by":self.name},layer=newlayer)
            else:
                # PyRat.Data.setAnnotation({"Created by":type(self)},layer=newlayer)
                PyRat.Data.setAnnotation({"Created by":type(self).__name__},layer=newlayer)

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
