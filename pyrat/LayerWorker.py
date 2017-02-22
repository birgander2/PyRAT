import logging
import pyrat


class LayerWorker(object):
    def __init__(self, shape=None, dtype=None, track=None, meta=None):

        self.attrs = {}
        self.attrs['_type'] = 'unknown'
        self.attrs['_shape'] = (0, 0)
        self.attrs['_lshape'] = (0,)
        self.attrs['_dshape'] = (0, 0)
        self.attrs['_dtype'] = 'unknown'
        self.attrs['_offset'] = (0, 0)
        self.data = None
        self.track = None
        self.prev = None
        if meta is not None:
            self.attrs.update(meta)
        self.name = pyrat.data.registerLayer(self)

    def run(self):
        pyrat.data.activateLayer(self.name)
        return self.name

    def setCrop(self, block, reset=False):
        if block[1] == 0:
            block[1] = self.attrs['_shape'][-2]
        if block[3] == 0:
            block[3] = self.attrs['_shape'][-1]
        if reset is True:
            self.attrs['_offset'] = (0, 0)
            self.attrs['_dshape'] = self.attrs['_shape'][-2:]
        else:
            self.attrs['_offset'] = (block[0], block[2])
            self.attrs['_dshape'] = (block[1] - block[0], block[3] - block[2])

    def deshape(self, shape):
        """
        Extracts layer shape and data shape from a numpy ndarray shape
        :returns: lshape, dshape
        """
        lshape = (1,)
        dshape = shape
        if len(shape) == 2:  # normal data
            pass
        elif len(shape) == 3:  # vector data
            lshape = (shape[0],)  # layer shape
            dshape = shape[1:]  # data shape
        elif len(shape) == 4:  # matrix data
            lshape = (shape[0], shape[1])  # layer shape
            dshape = shape[2:]  # data shape
        else:
            logging.error('Something wrong with array dimensions!')
            lshape = False
            dshape = False
        return lshape, dshape

    def setData(self):
        logging.error('Layer does not allow writing data')

    def getData(self):
        logging.error('Layer does not allow reading data')

    def setMeta(self):
        logging.error('Layer does not allow writing meta information')

    def getMeta(self):
        logging.error('Layer does not allow reading meta information')

    def setTrack(self):
        logging.error('Layer does not allow writing track data')

    def getTrack(self):
        logging.error('Layer does not allow reading track data')

    def info(self):
        print(self.attrs)
    
