import logging
import PyRat
import numpy as np

class LayerWorker(object):
    def __init__(self, shape=None, dtype=None, track=None, meta=None):
        
        self.attrs = {}
        self.attrs['_shape']  = (0,0)
        self.attrs['_lshape'] = (0,)
        self.attrs['_dshape'] = (0,0)
        self.attrs['_dtype']  = None
        self.data  = None
        self.track = None
        self.prev  = None
        if meta != None:
            self.attrs.update(meta)
        self.name = PyRat.Data.registerLayer(self)
     
    def run(self):
        PyRat.Data.activateLayer(self.name)
        return self.name
        
    def deshape(self, shape):
        """
        Extracts layer shape and data shape from a numpy ndarray shape
        :returns: lshape, dshape
        """
        lshape = (1,)
        dshape = shape
        if len(shape) == 2:                                               # normal data
            pass
        elif len(shape) == 3:                                             # vector data
            lshape = (shape[0],)                                          # layer shape
            dshape = shape[1:]                                            # data shape
        elif len(shape) == 4:                                             # matrix data
            lshape = (shape[0],shape[1])                                  # layer shape
            dshape = shape[2:]                                            # data shape
        else:
            logging.error('Something wrong with array dimensions!')
            lshape = False
            dshape = False
        return lshape, dshape
    
    def setData(self):
        logging.debug('Layer does not allow writing data')
        logging.error('Layer name unknown') 
    def getData(self):
        logging.debug('Layer does not allow reading data')
    
    def setMeta(self):
        logging.debug('Layer does not allow writing meta information')
    
    def getMeta(self):
        logging.debug('Layer does not allow reading meta information')
    
    def setTrack(self):
        logging.debug('Layer does not allow writing track data')
    
    def getTrack(self):
        logging.debug('Layer does not allow reading track data')
    
    def info(self):
        print self.attrs
    
 