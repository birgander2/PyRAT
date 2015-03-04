from __future__ import print_function
import h5py, os, logging, tempfile
import numpy as np
from pyrat.tools import deshape


class LayerData():
    def __init__(self, dir, *args):
        self.tmpdir = dir
        self.laynam = 1  # actual (unique) internal layer name
        self.layers = {}  # internal group / layer names
        self.active = []  # list of active layers
        self.shape = None  # data layer shape
        self.dtype = None  # data layer dtype
        self.lshape = None  # layer shape, e.g. (4,4) or (3,)
        self.dshape = None  # 2D data set shape, e.g. (2000,3000)
        self.ndim = None  # of dimensions, 2=image, 3=vector data, 4= matrix data

    def addLayer(self, array=None, shape=None, dtype='float32', memory=False, block='D', **kwargs):
        """
        Adds a new layer to the existing PyRat Data object. One can
        add one or multiple existing ndarrays (as list), or specify the 
        shape / dtype of a (single) new layer to be later filled with data.
        
        :author: Andreas Reigber
        :returns: list of layer names as strings
        """
        if isinstance(array, np.ndarray) or array is None:
            array = [array]
        if isinstance(array, tuple):
            array = list(array)

        addedgroups = []
        for arr in array:  # Loop over all arrays in list
            if arr is not None:
                shape = np.squeeze(arr).shape
                dtype = arr.dtype
            group = '/L' + str(self.laynam)
            if memory is False:
                logging.debug('Creating disc layer ' + group + ' = ' + str(dtype) + ' ' + str(shape))
                filename = tempfile.mktemp(suffix='.hd5', prefix='pyrat_', dir=self.tmpdir)
                self.layers[group] = DiscLayer(filename, group, shape, dtype, block=block)
            else:
                logging.debug('Creating memory layer ' + group + ' = ' + str(dtype) + ' ' + str(shape))
                self.layers[group] = MemoryLayer(group, shape, dtype, block=block)
            if arr is not None:
                lshape, dshape = deshape(shape)
                nchannels = np.prod(lshape)
                arr = arr.reshape((nchannels,) + dshape)
                self.layers[group].setData(arr)  # write data
            addedgroups.append(group)
            self.laynam += 1

        if len(addedgroups) == 1:
            return addedgroups[0]
        else:
            return addedgroups

    def activateLayer(self, layers, silent=False):
        """
        Activates a given layer or a list of layers or datasets. 
        """
        valid = self.existLayer(layers)
        if not isinstance(layers, list):
            layers = [layers, ]
            valid = [valid, ]
        elif len(layers) == 1:
            valid = [valid, ]
        self.active = [layer for (layer, val) in zip(layers, valid) if val]
        if len(self.active) == 1:
            self.active = self.active[0]

        if silent is False:
            logging.info('Activating ' + str(self.active))
        else:
            logging.debug('Activating ' + str(self.active))

        dtype = []
        ndim = []
        dshape = []
        lshape = []
        shape = []
        offset = []
        for layer in layers:
            group = '/' + layer.split('/')[1]
            s_dtype = self.layers[group].attrs["_dtype"]
            s_lshape = self.layers[group].attrs["_lshape"]
            s_dshape = self.layers[group].attrs["_dshape"]
            s_shape = self.layers[group].attrs["_shape"]
            s_offset = self.layers[group].attrs["_offset"]
            if 'D' in layer:
                s_lshape = (1, 0)
                s_shape = s_dshape
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
            offset.append(s_offset)

        if len(layers) == 1:
            self.dtype = dtype[0]
            self.dshape = dshape[0]
            self.lshape = lshape[0]
            self.shape = shape[0]
            self.ndim = ndim[0]
            self.offset = offset[0]
        else:
            self.dtype = dtype
            self.dshape = dshape
            self.lshape = lshape
            self.shape = shape
            self.ndim = ndim
            self.offset = offset

    def queryLayer(self, layers):
        """
        Returns a dict with layer information
        """
        valid = self.existLayer(layers)
        if not isinstance(layers, list):
            layers = [layers, ]
            valid = [valid, ]
        elif len(layers) == 1:
            valid = [valid, ]
        active = [layer for (layer, val) in zip(layers, valid) if val]
        if len(active) == 1:
            active = active[0]

        query = []
        for layer in layers:
            group = '/' + layer.split('/')[1]
            s_dtype = self.layers[group].attrs["_dtype"]
            s_lshape = self.layers[group].attrs["_lshape"]
            s_dshape = self.layers[group].attrs["_dshape"]
            s_shape = self.layers[group].attrs["_shape"]
            s_offset = self.layers[group].attrs["_offset"]
            if 'D' in layer:
                s_lshape = (1, 0)
                s_shape = s_dshape
            s_ndim = 3
            if s_lshape[0] == 1:
                s_ndim = 2
            elif len(s_lshape) == 2:
                s_ndim = 4
            query.append({'dtype': s_dtype, 'ndim': s_ndim, 'dshape': s_dshape, 'lshape': s_lshape, 'shape': s_shape})

        if len(layers) == 1:
            return query[0]
        else:
            return query

    def setData(self, array, layer=None, block=(0, 0, 0, 0)):
        """
        Fills a data layer with data. Gets called by addLayer method.
        """
        if layer is None:
            layer = self.active
        if isinstance(layer, str):
            layers = [layer, ]
            blocks = [block, ]
            array = [array, ]
        else:
            layers = layer
            if isinstance(block, tuple):
                blocks = [block] * len(layer)
            else:
                blocks = block
            if isinstance(array, np.ndarray):
                array = [array, ]

        for arr, layer, bl in zip(array, layers, blocks):
            valid = self.existLayer(layer)
            if valid is True:
                group = '/' + layer.split('/')[1]
                dshape = tuple(self.layers[group].attrs["_dshape"])
                dtype = self.layers[group].attrs['_dtype']
                if dtype != arr.dtype:
                    logging.error('dtype not compatible for layer ' + layer)
                    return
                if block[1] == 0: block = (block[0], self.layers[group].attrs['_dshape'][0], block[2], block[3])
                if block[3] == 0: block = (block[0], block[1], block[2], self.layers[group].attrs['_dshape'][1])
                bshape = (block[1] - block[0], block[3] - block[2], )
                if bshape != arr.shape[-2:]:
                    logging.error('array dimensions not compatible for layer ' + layer)
                    return
                self.layers[group].setData(arr, block=block, layer=layer)

    def getData(self, block=(0, 0, 0, 0), layer=None):
        if layer is None:
            layer = self.active
        if isinstance(layer, str):
            layers = [layer]
            blocks = [block]
        else:
            layers = layer
            if isinstance(block, tuple):
                blocks = [block] * len(layer)
            else:
                blocks = block
        array = []
        for layer, bl in zip(layers, blocks):
            valid = self.existLayer(layer)
            if valid is True:
                logging.debug('Reading layer ' + layer)
                group = '/' + layer.split('/')[1]
                dshape = tuple(self.layers[group].attrs["_dshape"])
                block = list(bl)
                array.append(self.layers[group].getData(block, layer=layer))
            else:
                array.append(np.zeros((0, 0)))
        if len(array) == 1:
            return array[0]
        else:
            return array

    def setAnnotation(self, annotation, layer=None):
        """
        Sets annotations (dict) to a layer or dataset.
        """
        if layer is None:
            layer = self.active
        if not isinstance(layer, list):
            layers = [layer]
        else:
            layers = layer
        if not isinstance(annotation, list):
            annotation = [annotation] * len(layers)
        for layer, anno in zip(layers, annotation):
            valid = self.existLayer(layer)
            if valid is True:
                group = '/' + layer.split('/')[1]
                self.layers[group].setMeta(anno, layer=layer)

    def getAnnotation(self, layer=None, key=False):
        """
        Returns all annotations of a layer or data set. 
        """
        if layer is None:
            layer = self.active
        layers = [layer] if not isinstance(layer, list) else layer
        annotation = []
        for layer in layers:
            valid = self.existLayer(layer)
            if valid is True:
                group = '/' + layer.split('/')[1]
                anno = self.layers[group].getMeta(layer=layer)
                if key is not False:
                    anno = anno[key]
                annotation.append(anno)
                # else:
                #annotation.append({})
        if len(annotation) == 1:
            return annotation[0]
        else:
            return annotation

    def calcBlock(self, block, layer=None):
        if layer is None:
            layer = self.active
        if isinstance(layer, list):
            layer = layer[0]
        group = '/' + layer.split('/')[1]
        offset = self.layers[group].attrs['_offset']
        dshape = self.layers[group].attrs['_dshape']

        block = list(block)
        block[0] += offset[0]
        block[2] += offset[1]
        if block[1] == 0:
            block[1] = dshape[0] + offset[0]
        else:
            block[1] += offset[0]
        if block[3] == 0:
            block[3] = dshape[1] + offset[1]
        else:
            block[3] += offset[1]
        return tuple(block)

    def resetCrop(self, layer=None):
        self.setCrop(reset=True, layers=layer)

    def setCrop(self, block=(0, 0, 0, 0), reset=False, layer=None):
        if layer is None:
            layer = self.active
        if isinstance(layer, str):
            layers = [layer, ]
            blocks = [block, ]
        else:
            layers = layer
            if isinstance(block, tuple):
                blocks = [block] * len(layer)
            else:
                blocks = block

        for layer, bl in zip(layers, blocks):
            valid = self.existLayer(layer)
            if valid is True:
                logging.info('Set crop to layer ' + layer)
                group = '/' + layer.split('/')[1]
                dshape = self.layers[group].attrs["_shape"][-2:]
                if block[0] < 0: block = (0, block[1], block[2], block[3])
                if block[1] > dshape[0]: block = (block[0], 0, block[2], block[3])
                if block[2] < 0: block = (block[0], block[1], 0, block[3])
                if block[3] > dshape[1]: block = (block[0], block[1], block[2], 0)
                self.layers[group].setCrop(block, reset=reset)
        dshape = []
        for layer in self.active:  # BUGGY LINES BELOW
            dshape.append(self.layers[group].attrs["_dshape"])
            self.dshape = dshape[0]
        if len(layers) == 1:
            self.dshape = tuple(dshape)
        else:
            self.dshape = tuple(dshape)

    def getLayerNames(self):
        return self.layers.keys()

    def getDataLayerNames(self, layer=None):
        if layer is None:
            layer = self.active
        layers = layer if isinstance(layer, list) else [layer]

        names = []
        for layer in layers:
            valid = self.existLayer(layer)
            if valid is True:
                group = '/' + layer.split('/')[1]
                nchannel = np.prod(self.layers[group].attrs['_lshape'])
                names.append([group + '/D' + str(channel) for channel in range(nchannel)])
        if len(names) == 1:
            names = names[0]
        return names

    def existLayer(self, layer):
        """
        Checks if layers are existing in data object. Returns boolean result.
        """
        layers = [layer] if not isinstance(layer, list) else layer

        valid = []
        for layer in layers:
            val = True
            try:
                foo = layer.split('/')
                group = '/' + foo[1]
                if group not in self.layers:
                    val = False
                else:
                    if 'D' in layer:
                        channel = int(foo[2][1:])
                        nchannel = np.prod(self.layers[group].attrs['_lshape'])
                        if channel < 0 or channel >= nchannel:
                            val = False
            except:
                val = False
            if val is False:
                logging.warning("Layer " + layer + " not existing!")
            valid.append(val)
        if len(valid) == 1:
            return valid[0]
        else:
            return valid

    def info(self):
        for group in ['/L' + str(l) for l in sorted([int(k[2:]) for k in self.layers.keys()])]:
            active = ' *' if group in self.active else ''
            print((self.layers[group].name + active).ljust(9), self.layers[group].attrs['_type'].ljust(15),
                  self.layers[group].attrs['_block'], self.layers[group].attrs['_dtype'].ljust(10),
                  self.layers[group].attrs['_shape'])

    def delLayer(self, layer, silent=False):
        """
        Deletes an entire layer from the object. Note that the memory (i.e. disc space) is not freed
        automatically for disc layers. A call to PyRat.Data.repack() is necessary, but it requires 
        copying around the data of all disc layers...
        """
        layers = layer if isinstance(layer, list) else [layer]
        for layer in layers:
            valid = self.existLayer(layer)
            if valid is True:
                if 'D' in layer:
                    logging.info('Cannot delete parts of layers')
                else:
                    if silent is False:
                        logging.info('Deletinging layer ' + layer)
                    else:
                        logging.debug('Deletinging layer ' + layer)

                    if self.layers[layer].attrs['_type'] == 'Disc':
                        self.layers[layer].group.close()
                        del self.layers[layer].group
                        os.remove(self.layers[layer].fn)
                    if layer in self.active:
                        self.active.remove(layer)   # TODO: Check for active data layers
                    del self.layers[layer]


class DiscLayer():
    def __init__(self, filename, group, shape, dtype, block='D', *args, **kwargs):
        self.fn = filename
        self.name = group
        self.group = h5py.File(self.fn, 'a')
        self.attrs = {}
        self.attrs['_type'] = 'Disc'
        lshape, dshape = deshape(shape)
        self.attrs['_shape'] = shape
        self.attrs['_lshape'] = lshape
        self.attrs['_dshape'] = dshape
        self.attrs['_offset'] = (0, 0)
        self.attrs['_dtype'] = str(dtype)
        self.attrs['_block'] = block
        self.group.create_group("P")  # Preview subgroup * not yet there
        self.group.create_dataset("D", (np.prod(lshape),) + dshape, dtype=dtype)  # create 2D Data layer

    def setCrop(self, block, reset=False):
        block = list(block)
        if block[1] == 0: block[1] = self.attrs['_shape'][-2]
        if block[3] == 0: block[3] = self.attrs['_shape'][-1]
        if reset is True:
            self.attrs['_offset'] = (0, 0)
            self.attrs['_dshape'] = self.attrs['_shape'][-2:]
        else:
            self.attrs['_offset'] = (block[0], block[2])
            self.attrs['_dshape'] = (block[1] - block[0], block[3] - block[2])
            # self.setMeta({'offset':self.attrs['_offset']})

    def setMeta(self, meta, layer=None):
        if meta is not None:
            for k, v in meta.items():
                if layer is not None and 'D' in layer and 'CH_' in k:
                    channel = int(layer.split('/')[2][1:])
                    if k not in self.group.attrs:
                        self.group.attrs[k] = [0] * np.prod(self.attrs['_lshape'])
                    ch_meta = list(self.group.attrs[k])
                    ch_meta[channel] = v
                    self.group.attrs[k] = ch_meta
                else:
                    if 'CH_' in k:
                        v = np.array(v)   # workaround: h5py doesn't support (yet?) unicode arrays
                        if 'U' in str(v.dtype):  # if unicode
                            v = v.astype('S')    # then convert to string
                    self.group.attrs[k] = v

    def getMeta(self, key=False, layer=None):
        meta = dict(self.group.attrs)
        for k, v in meta.items():
            if isinstance(v, np.ndarray) and (v.dtype in ['|S2', '|S1', '|S5']):
                meta[k] = [foo.decode() for foo in v]
        if layer is not None and 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            ch_meta = [k for k in meta.keys() if 'CH_' in k]
            for key in ch_meta:
                meta[key] = meta[key][channel]
        for k in meta.keys():
            if k[0] == '_':
                del meta[k]
        return meta

    def setData(self, array, block=(0, 0, 0, 0), layer=None):
        offset = self.attrs['_offset']
        block = list(block)
        block[0] += offset[0]
        block[2] += offset[1]
        if block[1] == 0:
            block[1] = self.attrs['_dshape'][0] + offset[0]
        else:
            block[1] += offset[0]
        if block[3] == 0:
            block[3] = self.attrs['_dshape'][1] + offset[1]
        else:
            block[3] += offset[1]

        if layer is None or layer == self.name:
            nchannels = np.prod(self.attrs["_lshape"])
            if self.attrs['_block'] == 'D':
                self.group["/D"][..., block[0]:block[1], block[2]:block[3]] = array.reshape(
                    (nchannels,) + array.shape[-2:])
            elif self.attrs['_block'] == 'T':
                self.group["/D"][..., block[0]:block[1], :] = array.reshape((nchannels,) + array.shape[-2:])
            elif self.attrs['_block'] == 'O':
                self.group["/D"][...] = array.reshape((nchannels,) + array.shape[-2:])
        elif 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            if self.attrs['_block'] == 'D':
                self.group["/D"][channel, block[0]:block[1], block[2]:block[3]] = array
            elif self.attrs['_block'] == 'T':
                self.group["/D"][channel, block[0]:block[1], :] = array
            elif self.attrs['_block'] == 'O':
                self.group["/D"][channel, ...] = array
        else:
            logging.error('Layer name unknown')
            stop()

    def getData(self, block=(0, 0, 0, 0), layer=None):
        offset = self.attrs['_offset']
        block = list(block)
        block[0] += offset[0]
        block[2] += offset[1]
        if block[1] == 0:
            block[1] = self.attrs['_dshape'][0] + offset[0]
        else:
            block[1] += offset[0]
        if block[3] == 0:
            block[3] = self.attrs['_dshape'][1] + offset[1]
        else:
            block[3] += offset[1]
        if layer is None or layer == self.name:
            lshape = tuple(self.attrs["_lshape"])
            dshape = tuple(self.attrs["_dshape"])
            if self.attrs['_block'] == 'D':
                bshape = (block[1] - block[0], block[3] - block[2])
                # print(block,lshape+bshape)
                return np.squeeze(
                    np.reshape(self.group["/D"][..., block[0]:block[1], block[2]:block[3]], lshape + bshape))
            elif self.attrs['_block'] == 'T':
                bshape = (block[1] - block[0], dshape[-1])
                return np.squeeze(np.reshape(self.group["/D"][..., block[0]:block[1], :], lshape + bshape))
            elif self.attrs['_block'] == 'O':
                bshape = tuple(self.attrs["_dshape"])
                return np.squeeze(np.reshape(self.group["/D"][...], lshape + bshape))
        elif 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            if self.attrs['_block'] == 'D':
                return np.squeeze(self.group["/D"][channel, block[0]:block[1], block[2]:block[3]])
            elif self.attrs['_block'] == 'T':
                return np.squeeze(self.group["/D"][channel, block[0]:block[1], :])
            elif self.attrs['_block'] == 'O':
                return np.squeeze(self.group["/D"][channel, ...])
        else:
            logging.error('Layer name unknown')
            stop()


class MemoryLayer():
    def __init__(self, name, shape, dtype, *args, **kwargs):
        self.name = name
        self.attrs = {}
        self.attrs['_type'] = 'Memory'
        lshape, dshape = deshape(shape)
        self.attrs['_shape'] = shape
        self.attrs['_lshape'] = lshape
        self.attrs['_dshape'] = dshape
        self.attrs['_offset'] = (0, 0)
        self.attrs['_dtype'] = str(dtype)
        self.attrs['_block'] = kwargs['block']
        nchannels = np.prod(lshape)
        self.data = np.empty((nchannels,) + dshape, dtype=dtype)

    def setCrop(self, block, reset=False):
        if block[1] == 0: block[1] = self.attrs['_shape'][-2]
        if block[3] == 0: block[3] = self.attrs['_shape'][-1]
        if reset is True:
            self.attrs['_offset'] = (0, 0)
            self.attrs['_dshape'] = self.attrs['_shape'][-2:]
        else:
            self.attrs['_offset'] = (block[0], block[2])
            self.attrs['_dshape'] = (block[1] - block[0], block[3] - block[2])

    def setMeta(self, meta, layer=None):
        for k, v in meta.items():
            if layer is not None and 'D' in layer and 'CH_' in k:
                channel = int(layer.split('/')[2][1:])
                if k not in self.attrs:
                    self.attrs[k] = [0] * np.prod(self.attrs['_lshape'])
                self.attrs[k][channel] = v
            else:
                self.attrs[k] = v

    def getMeta(self, key=False, layer=None):
        meta = self.attrs.copy()
        if 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            ch_meta = [k for k in meta.keys() if 'CH_' in k]
            for key in ch_meta:
                meta[key] = meta[key][channel]
        for k in meta.keys():
            if k[0] == '_':
                del meta[k]
        return meta

    def setData(self, array, block=(0, 0, 0, 0), layer=None):
        if self.attrs['_block'] is False: block = (0, 0, 0, 0)
        offset = self.attrs['_offset']
        block = list(block)
        block[0] += offset[0]
        block[2] += offset[1]
        if block[1] == 0:
            block[1] = self.attrs['_dshape'][0] + offset[0]
        else:
            block[1] += offset[0]
        if block[3] == 0:
            block[3] = self.attrs['_dshape'][1] + offset[1]
        else:
            block[3] += offset[1]

        if layer is None or layer == self.name:
            nchannels = np.prod(self.attrs["_lshape"])
            self.data[..., block[0]:block[1], block[2]:block[3]] = array.reshape((nchannels,) + array.shape[-2:])
        elif 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            self.data[channel, block[0]:block[1], block[2]:block[3]]
        else:
            logging.error('Layer name unknown')

    def getData(self, block=(0, 0, 0, 0), layer=None):
        if self.attrs['_block'] is False: block = (0, 0, 0, 0)
        offset = self.attrs['_offset']
        block = list(block)
        block[0] += offset[0]
        block[2] += offset[1]
        if block[1] == 0:
            block[1] = self.attrs['_dshape'][0] + offset[0]
        else:
            block[1] += offset[0]
        if block[3] == 0:
            block[3] = self.attrs['_dshape'][1] + offset[1]
        else:
            block[3] += offset[1]

        if layer is None or layer == self.name:
            bshape = (block[1] - block[0], block[3] - block[2])
            lshape = tuple(self.attrs["_lshape"])
            return np.squeeze(np.reshape(self.data[..., block[0]:block[1], block[2]:block[3]], lshape + bshape))
        elif 'D' in layer:
            channel = int(layer.split('/')[2][1:])
            return np.squeeze(self.data[channel, block[0]:block[1], block[2]:block[3]])
        else:
            logging.error('Layer name unknown') 

