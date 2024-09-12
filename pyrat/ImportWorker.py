import pyrat
import logging
from .tools import AttrDict


class ImportWorker(pyrat.Worker):
    def __init__(self, *args, **kwargs):
        super(ImportWorker, self).__init__(*args, **kwargs)
        self.nthreads = 1
        self.block = 'D'

    def main(self, *args, **kwargs):
        # Enabling import of multiple files at once
        if hasattr(self, 'file'):
            if type(self.file) not in (list, tuple):
                files = [self.file]
            else:
                files = self.file
        else:
            files = [None]

        layers = list()
        for self.file in files:
            size = tuple(self.getsize(*args, **kwargs))
            size += (1,)*(2-len(size))
            if size[0] is None:                  # getsize not overloaded -> use full image import
                data, meta = self.reader(*args, **kwargs)
                if data is None and meta is None:
                    logging.warning("Nothing imported!!!")
                    return False
                if data is None:
                    logging.debug("No image data imported")
                if meta is None:
                    logging.debug("No meta data imported")
                if isinstance(meta, dict):        # AttrDict is preferred over normal dict
                    meta = AttrDict(meta)
                newlayer = pyrat.data.addLayer(array=data, block=self.block)
                pyrat.data.setAnnotation(meta, layer=newlayer)
                pyrat.data.activateLayer(newlayer)
                layers.append(newlayer)
            elif size[0] is False:                # getsize failed in some sense -> return False
                return False
            else:                                 # blockwise import
                newlayer = self.layer_fromfunc(self.block_reader, size=size, silent=False)
                meta = self.getmeta(*args, **kwargs)
                if isinstance(meta, dict):        # AttrDict is preferred over normal dict
                    meta = AttrDict(meta)
                self.close(*args, **kwargs)
                if meta is not False:
                    pyrat.data.setAnnotation(meta, layer=newlayer)
                else:
                    logging.debug("No meta data imported")
                pyrat.data.activateLayer(newlayer)
                layers.append(newlayer)
        if len(layers) == 1:
            return layers[0]
        else:
            return layers

    def getmeta(self, *args, **kwargs):
        return False

    def getsize(self, *args, **kwargs):
        return None, None

    def reader(self, *args, **kwargs):
        logging.warning("No image data imported")
        return False

    def block_reader(self, *args, **kwargs):
        return False

    def close(self, *args, **kwargs):
        """
        Close the file here (to be overloaded)
        """
        pass

