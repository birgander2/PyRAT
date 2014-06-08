import h5py, os, logging, pdb
import numpy as np
import itertools
#import line_profiler

class LayerData(h5py.File):
    def __init__(self, filename, *args):
        read_flag = False
        if os.path.exists(filename):
            read_flag = True
        super(LayerData, self).__init__(filename, 'a', *args)
        self.laynam = 1                           # actual (unique) internal layer name
        self.layers = {}                          # internal group / layer names
        self.active = []                          # list of active layers
        self.shape  = None                        # data layer shape
        self.dtype  = None                        # data layer dtype
        self.lshape = None                        # layer shape, e.g. (4,4) or (3,)
        self.dshape = None                        # 2D data set shape, e.g. (2000,3000)
        self.ndim   = None                        # of dimensions, 2=image, 3=vector data, 4= matrix data

    def registerLayer(self, layer, name=None):
        if name == None:
            name = '/L'+str(self.laynam)
            self.laynam += 1
        self.layers[name] = layer
        logging.info('Registering layer '+name)
        return name

    def addLayer(self, array=None, shape=None, dtype='float32', replace=False, memory=False):
        """
        Adds a new layer to the existing PyRat Data object. One can
        add one or multiple existing ndarrays (as list), or specify the
        shape / dtype of a (single) new layer to be later filled with data.

        :author: Andreas Reigber
        :returns: list of layer names as strings
        """
        if isinstance(array,np.ndarray) or array == None:
            array = (array,)

        if replace == True:
            layers = self.active
            if not isinstance(layers,tuple):
                layers = (layers,)
            for arr, layer in zip(array, layers):
                if arr != None:
                    shape = arr.shape
                    dtype = arr.dtype
                lshape, dshape = deshape(shape)
                nchannels = np.prod(lshape)
                if arr != None:
                    arr = arr.reshape((nchannels,)+dshape)
                ds = layer
                logging.debug('Updating content of layer '+ds)
                if arr != None:
                    self.setData(arr, layer=str(ds))
            return self.active

        addedgroups = []
        for arr in array:                                               # Loop over all arrays in list
            if arr != None:
                shape = arr.shape
                dtype = arr.dtype
            group = '/L'+str(self.laynam)
            if memory == False:
                logging.info('Creating disc layer '+group+' = '+str(dtype)+' '+str(shape))
                grp   = self.create_group(group)
                self.layers[group] = DiscLayer(grp, shape, dtype)
            else:
                logging.info('Creating memory layer '+group+' = '+str(dtype)+' '+str(shape))
                self.layers[group] = MemoryLayer(group, shape, dtype)
            if arr != None:
                lshape, dshape = deshape(shape)
                nchannels = np.prod(lshape)
                arr = arr.reshape((nchannels,)+dshape)
                self.layers[group].setData(arr)                               # write data
            addedgroups.append(group)
            self.laynam += 1
        addedgroups = tuple(addedgroups)
        if len(addedgroups) == 1:
            return addedgroups[0]
        else:
            return addedgroups

    def setData(self, array, layer=None, block=None):
        """
        Fills a data layer with data. Gets called by addLayer method.
        """
        if layer == None:
            layer = self.active
        if block == None:
            block = [0,0,0,0]
        if isinstance(layer,str):
            layers = (layer,)
            blocks = (block,)
            array  = (array,)
        else:
            layers = layer
            if isinstance(block,list):
                blocks = tuple([block]*len(layer))
            else:
                blocks = block

        for arr, layer, bl in zip(array, layers, blocks):
            group = '/'+layer.split('/')[1]
            if group not in self.layers:
                logging.error('Layer '+group+' not existing')
                pdb.set_trace()
            dshape = tuple(self.layers[group].attrs["_dshape"])

            block = list(bl)
            self.layers[group].setData(arr, block=block, layer=layer)

    def activateLayer(self, layers):
        """
        Activates a given layer or a list of layers or datasets.
        """

        if isinstance(layers,tuple):
            if not all(['/'+layer.split('/')[1] in self.layers for layer in layers]):
                logging.error("At least one layer is not existing!")
                pdb.set_trace()
            #if len(layers) == 1:
                #layers = layers[0]
        else:
            if not '/'+layers.split('/')[1] in self.layers:
                logging.error("Layer is not existing!")
                pdb.set_trace()

        self.active = layers
        logging.info('Activating '+str(layers))

        if not isinstance(layers,tuple):
            layers = (layers,)

        dtype  = []
        ndim   = []
        dshape = []
        lshape = []
        shape  = []
        for layer in layers:
            group = '/'+layer.split('/')[1]
            s_dtype  = self.layers[group].attrs["_dtype"]
            s_lshape = self.layers[group].attrs["_lshape"]
            s_dshape = self.layers[group].attrs["_dshape"]
            s_shape  = self.layers[group].attrs["_shape"]
            if 'D' in layer:
                s_lshape = (1,0)
                s_shape  = s_dshape
            s_ndim = 3
            if s_lshape[0] == 1:
                s_ndim = 2
            elif len(s_lshape) == 2:
                s_ndim = 4
            dtype.append(s_dtype)
            ndim.append(s_ndim)
            dshape.append(s_dshape)
            lshape.append(s_lshape)
            shape.append(s_shape)

        if len(layers) == 1:
            self.dtype  = dtype[0]
            self.dshape = dshape[0]
            self.lshape = lshape[0]
            self.shape  = shape[0]
            self.ndim   = ndim[0]
        else:
            self.dtype  = tuple(dtype)
            self.dshape = tuple(dshape)
            self.lshape = tuple(lshape)
            self.shape  = tuple(shape)
            self.ndim   = tuple(ndim)


    def getData(self, block=None, layer=None):
        if layer == None:
            layer = self.active
        if block == None:
            block = [0,0,0,0]
        if isinstance(layer,str):
            layers = (layer,)
            blocks = (block,)
        else:
            layers = layer
            if isinstance(block,list):
                blocks = tuple([block]*len(layer))
            else:
                blocks = block

        array = []
        for layer, bl in zip(layers, blocks):
            logging.debug('Reading layer '+layer)
            group = '/'+layer.split('/')[1]
            dshape = tuple(self.layers[group].attrs["_dshape"])
            block = list(bl)
            array.append(self.layers[group].getData(block, layer=layer))

        if len(array) == 1:
            return array[0]
        else:
            return tuple(array)

    def delLayer(self, layers):
        """
        Deletes an entire layer from the object. Note that the memory (i.e. disc space) is not freed
        automatically. A call to PyRat.cleanup() is necessary to do so, but this requires copying around
        all the data of the object.
        """
        if not isinstance(layers,tuple):
            layers = (layers,)
        for layer in layers:
            if layer in self.layers:
                logging.info('Deletinging layer '+layer)
                del self.layers[layer]
            else:
                logging.error('Layer '+str(layer)+' not existing')
                pdb.set_trace()

    def setAnnotation(self, annotation, layers=False):
        """
        Sets annotations (dict) to a layer or dataset.
        """
        if layers == False:
            layers = self.active
        if not isinstance(layers,tuple):
            layers = (layers,)
        if not isinstance(annotation, tuple):
            annotation = tuple([annotation]*len(layers))
        for layer, anno in zip(layers, annotation):
            if layer in self.layers:
                self.layers[layer].setMeta(anno)
            else:
                logging.error('Layer '+str(layer)+' not existing')
                pdb.set_trace()

    def getAnnotation(self, layers=False, key=False):
        """
        Returns all annotations of a layer or data set.
        """
        if layers == False:
            layers = self.active
        if not isinstance(layers,tuple):
            layers = (layers,)
        annotation = []
        for layer in layers:
            if layer in self.layers:
                anno = self.layers[layer].getMeta()
                #pdb.set_trace()
                if key != False:
                    anno = anno[key]
                annotation.append(anno)
            else:
                logging.error('Layer '+str(layer)+' not existing')
                pdb.set_trace()
                annotation.append({})
        if len(annotation) == 1:
            return annotation[0]
        else:
            return tuple(annotation)

    def setBlock(self, block=[0,0,0,0]):
        if isinstance(self.dshape, tuple):
            logging.error('Multiple layers selected!')
            self.block = block

        if block[0] < 0: block[0] = 0
        if block[0] >= self.dshape[0]: block[0] = self.dshape[0]-1
        if block[1] < 0: block[1] = 0
        if block[1] > self.dshape[0]: block[1] = self.dshape[0]
        if block[2] < 0: block[2] = 0
        if block[2] >= self.dshape[1]: block[2] = self.dshape[1]-1
        if block[3] < 0: block[3] = 0
        if block[3] > self.dshape[1]: block[3] = self.dshape[1]
        logging.info('Data block selected: y='+str(block[0])+'->'+str(block[1])+', x='+str(block[2])+'->'+str(block[3]))
        self.block = block


    def setTrack(self, track, layer=False):
        #logging.debug('setTrack not implemented')
        pass

    def getTrack(self, layer=False):
        return None
        #logging.debug('getTrack not implemented')

    def info(self):
        def printname(name):
            print name
        self.visit(printname)
        print self.layers

class DiscLayer():
    def __init__(self, hdfgroup, shape, dtype, *args, **kwargs):
        self.group = hdfgroup
        self.name  = hdfgroup.name
        self.attrs = {}
        self.attrs['_type']   = 'HDF5'
        lshape, dshape = deshape(shape)
        self.attrs['_shape']  = shape
        self.attrs['_lshape'] = lshape
        self.attrs['_dshape'] = dshape
        self.attrs['_dtype']  = str(dtype)
        self.group.create_group("P")                                          # Preview subgroup * not yet there
        self.group.create_dataset("D", (np.prod(lshape),)+dshape, dtype=dtype)      # create 2D Data layer
        self.group.create_dataset("T",(dshape[0], 4 ), dtype='float64')       # track data

    def setMeta(self, meta):
        for k,v in meta.items():
            self.group.attrs[k] = v

    def getMeta(self, key=False):
        meta = dict(self.group.attrs)
        for k,v in meta.items():
            if isinstance(v,np.ndarray) and v.dtype == '|S2': meta[k] = list(v)
        for k in meta.keys():
            if k[0] == '_':
                del meta[k]
        return meta

    def setData(self, array, block=[0,0,0,0], layer=None):
        if block[1] == 0: block[1] =  self.attrs['_dshape'][0]
        if block[3] == 0: block[3] =  self.attrs['_dshape'][1]
        if layer == None or layer == self.name:
            nchannels = np.prod(self.attrs["_lshape"])
            self.group[self.name+"/D"][...,block[0]:block[1],block[2]:block[3]] = array.reshape((nchannels,)+array.shape[-2:])
        elif 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            self.group[self.name+"/D"][channel,block[0]:block[1],block[2]:block[3]]
        else:
            logging.error('Layer name unknown')
            pdb.set_trace()

    def getData(self, block=[0,0,0,0], layer=None):
        if block[1] == 0: block[1] =  self.attrs['_dshape'][0]
        if block[3] == 0: block[3] =  self.attrs['_dshape'][1]
        if layer == None or layer == self.name:
            bshape = (block[1]-block[0],block[3]-block[2])
            lshape = tuple(self.attrs["_lshape"])
            return np.squeeze(np.reshape(self.group[self.name+"/D"][...,block[0]:block[1],block[2]:block[3]],lshape+bshape))
        elif 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            return np.squeeze(self.group[self.name+"/D"][channel,block[0]:block[1],block[2]:block[3]])
        else:
            logging.error('Layer name unknown')
            pdb.set_trace()

class MemoryLayer():

    def __init__(self, name, shape, dtype, *args, **kwargs):
        self.name  = name
        self.attrs = {}
        self.attrs['_type']   = 'Numpy'
        lshape, dshape = deshape(shape)
        self.attrs['_shape']  = shape
        self.attrs['_lshape'] = lshape
        self.attrs['_dshape'] = dshape
        self.attrs['_dtype']  = str(dtype)
        nchannels = np.prod(lshape)
        self.data = np.empty((nchannels,)+dshape, dtype=dtype)

    def setMeta(self, meta):
        for k,v in meta.items():
            self.attrs[k] = v

    def getMeta(self, key=False):
        meta = self.attrs.copy()
        for k in meta.keys():
            if k[0] == '_':
                del meta[k]
        return meta

    def setData(self, array, block=[0,0,0,0], layer=None):
        if block[1] == 0: block[1] =  self.attrs['_dshape'][0]
        if block[3] == 0: block[3] =  self.attrs['_dshape'][1]
        if layer == None or layer == self.name:
            nchannels = np.prod(self.attrs["_lshape"])
            self.data[...,block[0]:block[1],block[2]:block[3]] = array.reshape((nchannels,)+array.shape[-2:])
        elif 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            self.data[channel,block[0]:block[1],block[2]:block[3]]
        else:
            logging.error('Layer name unknown')
            pdb.set_trace()

    def getData(self, block=[0,0,0,0], layer=None):
        if block[1] == 0: block[1] =  self.attrs['_dshape'][0]
        if block[3] == 0: block[3] =  self.attrs['_dshape'][1]
        if layer == None or layer == self.name:
            bshape = (block[1]-block[0],block[3]-block[2])
            lshape = tuple(self.attrs["_lshape"])
            return np.squeeze(np.reshape(self.data[...,block[0]:block[1],block[2]:block[3]],lshape+bshape))
        elif 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            return np.squeeze(self.data[channel,block[0]:block[1],block[2]:block[3]])
        else:
            logging.error('Layer name unknown')
            pdb.set_trace()

def deshape(shape):
    """
    Extracts layer shape and data shape from a numpy ndarray shape
    :returns: lshape, dshape
    """
    lshape = (1,)
    dshape = shape
    if len(shape) == 2:                                               # normal data
        pass
    elif len(shape) == 3:                                               # vector data
        lshape = (shape[0],)                                          # layer shape
        dshape = shape[1:]                                            # data shape
    elif len(shape) == 4:                                               # matrix data
        lshape = (shape[0],shape[1])                                  # layer shape
        dshape = shape[2:]                                            # data shape
    else:
        logging.error('Something wrong with array dimensions!')
        lshape = False
        dshape = False
    return lshape, dshape


