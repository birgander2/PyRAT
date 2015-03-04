from __future__ import print_function
import pyrat
import logging
from PyQt4 import QtGui, QtCore


class ImportWorker(pyrat.Worker):
    def __init__(self, *args, **kwargs):
        super(ImportWorker, self).__init__(*args, **kwargs)
        self.nthreads = 1
        self.block = 'D'

    def run(self, *args, **kwargs):

        para = [foo['var'] for foo in self.para]
        self.checkpara(kwargs, para)
        logging.info(
            self.name + '  ' + str(dict((k, v) for k, v in self.__dict__.items() if k in para or k in kwargs)))

        size = tuple(self.getsize(*args, **kwargs))
        if size[0] is None:                  # getsize not overloaded -> use full image import
            data, meta = self.reader(*args, **kwargs)
            if data is None and meta is None:
                logging.warning("Nothing imported!!!")
                return False
            if data is None:
                logging.warning("No image data imported")
            if meta is None:
                logging.warning("No meta data imported")
            newlayer = pyrat.data.addLayer(data, block=self.block, *args, **kwargs)
            pyrat.data.setAnnotation(meta, layer=newlayer)
            pyrat.data.activateLayer(newlayer)
            return newlayer
        elif size[0] is False:                # getsize failed in some sense -> return False
            return False
        else:                                 # blockwise import
            newlayer = self.layer_fromfunc(self.block_reader, size=size, silent=False)
            meta = self.getmeta(*args, **kwargs)
            if meta is not False:
                pyrat.data.setAnnotation(meta, layer=newlayer)
            else:
                logging.warning("No meta data imported")
            pyrat.data.activateLayer(newlayer)
            return newlayer

    def getmeta(self, *args, **kwargs):
        return False

    def getsize(self, *args, **kwargs):
        return None, None

    def reader(self, *args, **kwargs):
        logging.warning("No image data imported")
        return False

    def block_reader(self, *args, **kwargs):
        return False

