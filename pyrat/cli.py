import h5py
import logging
import tempfile
import numpy as np
import pickle
import pyrat
from .tools import AttrDict


def listlayer(*args, **kwargs):
    """
    Prints a list of existing layers. The active layers are marked with a star.
    """
    from pyrat import data

    data.listLayer()


info = listlayer


def active(*args, **kwargs):
    """
    Returns the names of the currently active layers.
    """
    from pyrat import data

    return data.active


def activate(*args, silent=False, **kwargs):
    """
    Activates one or several layers.
    :param layer: String with layer name to activate (alternatively: list of layer names)
    :param silent: stay silent, do not print output
    """
    from pyrat import data

    if 'layer' in kwargs:
        layer = kwargs['layer']
    elif len(args) > 0:
        layer = args[0]

    if len(layer) > 0:
        data.activateLayer(layer, silent=silent)


def adddata(array, block='D', **kwargs):
    """
    Generates a new pyrat layer from a numpy ndarray
    :param array: The numpy array to be saved in the new layer
    :param block: "D" for normal data (default), "T" for tracks
    """
    from pyrat import data

    layer = data.addLayer(array=array, block=block)
    activate(layer)
    return layer


addlayer = adddata


def updatelayer(array, **kwargs):
    """
    Updates the content of a pyrat layer with the content of a provided ndarray.
    Warning: Do not change dimensions or dtype!
    :param array: The numpy array to be written in the new layer
    :param layer: Optional layer name to use as target (default: active layer)
    """
    from pyrat import data

    data.setData(array, **kwargs)


def delete(*args, silent=False, **kwargs):
    """
    Deletes one or several layers.
    :param layer: String with layer name to delete (alternatively: list of layer names)
    :param silent: stay silent, do not print output
    """
    from pyrat import data

    if 'layer' in kwargs:
        layer = kwargs['layer']
    elif len(args) > 0:
        layer = args[0]
    else:
        layer = active()

    if len(layer) > 0:
        data.delLayer(layer, silent=silent)


def query(*args, **kwargs):
    """
    Query the content / data structure of a layer

    :param layer: String with layer name to delete (alternatively: list of layer names)
    """
    from pyrat import data

    if 'layer' in kwargs:
        pass
    elif len(args) > 0:
        kwargs['layer'] = args[0]
    else:
        kwargs['layer'] = active()

    if isinstance(kwargs['layer'], str) or len(kwargs['layer']) == 1:
        val = data.queryLayer(**kwargs)
        if val is None:
            return
        else:
            return AttrDict(val)
    else:
        logging.error("This works only for a single layer (up to now)")


def getdata(*args, **kwargs):
    """
    Copies the content of a layer into a numpy array
    :param layer: layer name to extract (optional, default=active)
    :return: layer content as numpy arry
    """
    from pyrat import data

    if 'layer' in kwargs:
        pass
    elif len(args) > 0:
        kwargs['layer'] = args[0]
    else:
        kwargs['layer'] = active()

    return data.getData(**kwargs)


def getmeta(*args, **kwargs):
    """
    Returns the meta data of a layer as dict
    :param layer: layer name to extract meta data (optional, default=active)
    :return: meta data
    """
    from pyrat import data

    if 'layer' in kwargs:
        pass
    elif len(args) > 0:
        kwargs['layer'] = args[0]
    else:
        kwargs['layer'] = active()

    return data.getAnnotation(**kwargs)


def setmeta(meta, **kwargs):
    """
    Adds metadata to a layer
    :param meta: meta data as dict
    :param layer: layer name to extract meta data (optional, default=active)
    """
    from pyrat import data

    data.setAnnotation(meta, **kwargs)


def setlayername(name, **kwargs):
    """
    Sets the name of a layer
    :param name: string with a descriptive name
    :param layer: layer identifier to assign the name (optional, default=active)
    """
    from pyrat import data

    data.setLayerName(name, **kwargs)


def expose(*args, **kwargs):
    """
    Expose the layer content for direct manipulation
    The returned object should behave as a numpy array (but might be of type hdf5). Manipulations of it
    are instantanously written to the layer. Do not change size or dtype of data!
    Attention: multidimensional arrays appear flattened!
    :param  layer name to extract meta data (optional, default=active)
    :return: layer object
    """
    from pyrat import data

    layer = []
    if 'layer' in kwargs:
        layer.append(kwargs['layer'])
    elif len(args) > 0:
        layer = args[0]
    else:
        layer = active()

    if len(layer) == 1:
        return data.exposeRaw(layer=layer)
    else:
        logging.error("This works only for a single layer (up to now)")


def freeze(filename, *args, **kwargs):
    """
    Freezes the entire pyrat session (the layer structure) for later usage.
    A list of local variables can be saved together with the layers.
    Example: freeze("session.hd5", var1, var2, var3)
    :param filename: File where to save all the data (might get very large!)
    :param args: List of local variables to be saved together with the layer
    """
    from pyrat import data

    file = h5py.File(filename, 'w')
    for layer in data.layers:
        if data.layers[layer].attrs['_type'] == 'Memory':
            print("WARNING: memory layers cannot be frozen...")
        else:
            inp = h5py.File(data.layers[layer].fn, 'r')
            print("Saving layer", layer, inp.filename)
            inp.flush()
            out = file.create_group(layer)
            data.layers[layer].file.close()
            del data.layers[layer].file
            del data.layers[layer].group
            pick = pickle.dumps(data.layers[layer])
            data.layers[layer].file = h5py.File(data.layers[layer].fn, 'a')
            data.layers[layer].group = data.layers[layer].file['/D']
            out.create_dataset("pickle", data=np.fromstring(pick, dtype=np.uint8))
            h5py.h5o.copy(inp.id, str.encode("D"), out.id, str.encode("D"))
    file.attrs["laynam"] = data.laynam
    if isinstance(data.active, list):
        file.attrs["active"] = np.array(data.active, dtype="S4")
    else:
        file.attrs["active"] = data.active
    file.create_dataset("locals", data=np.fromstring(pickle.dumps(args), dtype=np.uint8))
    file.close()


def defreeze(filename, *args, **kwargs):
    """
    Restores are previously session saved with 'freeze'
    Example:  var1, var2, var3 = defreeze("session.hd5")
    :param filename:  Session file saved with 'freeze'
    :return: Additional variables stored during the freeze process.
    """
    from pyrat import data

    file = h5py.File(filename, 'r')
    data.layers = {}
    local_var = ()
    for layer in file:
        if layer != 'locals':
            fn = tempfile.mktemp(suffix='.hd5', prefix='pyrat_', dir=data.tmpdir)
            print("Reading layer", layer, fn)
            out = h5py.File(fn, 'w')
            h5py.h5o.copy(file.id, str.encode(layer + '/D'), out.id, str.encode("D"))
            obj = pickle.loads((bytes(file[layer]['pickle'][...])))
            data.layers['/' + layer] = obj
            data.layers['/' + layer].fn = fn
            data.layers['/' + layer].file = h5py.File(fn, 'a')
            data.layers['/' + layer].group = data.layers['/' + layer].file['/D']
        else:
            local_var = pickle.loads((bytes(file['locals'][...])))
    data.laynam = file.attrs["laynam"]
    active = file.attrs["active"]
    if isinstance(active, str):
        data.active = active
    else:
        data.active = list(active.astype("U"))
    file.close()
    data.activateLayer(data.active, silent=True)

    if len(local_var) > 0:
        return local_var


def show(*args, method='amplitude', **kwargs):
    """
    Launches the PyRat GUI
    """

    import pyrat
    import sys
    from PyQt5 import QtWidgets, QtCore

    if method not in ['amplitude', 'intensity', 'phase', '0.0->1.0', 'min->max', 'labels']:
        method = 'amplitude'

    if 'layer' in kwargs:
        activate(kwargs['layer'])
    elif len(args) == 1:
        activate(args[0])

    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)
        new_app = True
    else:
        new_app = False
    pyrat.app = pyrat.viewer.MainWindow()
    pyrat.app.updateViewer(method=method)
    if 'palette' in kwargs:
        pyrat.app.config['palette'] = kwargs['palette']
    pyrat.app.updateViewer(method=method)
    if new_app is True:
        app.exec_()

gui = show

def help(*args, **kwargs):
    import sys, types
    from inspect import getmembers, isfunction, ismodule

    if len(args) > 0:
        for arg in args:
            if ismodule(arg) and hasattr(arg, 'help'):
                logging.info(arg.help.__doc__)
                arg.help()
            else:
                logging.info(arg.__doc__)

    else:
        current_module = sys.modules[__name__]
        modules = getmembers(current_module, isfunction)
        logging.info("")
        logging.info("Available main level functions:")
        logging.info("")
        for mod in modules:
            if 'pyrat' in mod[1].__module__ and \
                    mod[1].__name__ not in ['help']:
                doc = str(mod[1].__doc__)
                if doc != 'None':
                    doc = doc.split('\n')[1].lstrip()
                else:
                    doc = "-"
                logging.info((mod[0] + "()").ljust(20) + doc)

        logging.info("")
        logging.info("Available pyrat modules:")
        logging.info("")
        modules = getmembers(pyrat, ismodule)
        modules.append(('plugins', pyrat.plugins))
        modules = sorted(modules, key=lambda x: x[0])

        for mod in modules:
            if 'pyrat' in mod[1].__name__ and hasattr(mod[1], 'help') and \
                    mod[1].__name__ not in ['pyrat.cli', 'pyrat']:
                doc = str(mod[1].help.__doc__)
                if doc != 'None':
                    doc = doc.split('\n')[1].lstrip()
                else:
                    doc = "-"
                logging.info(mod[0].ljust(20) + doc)
        logging.info("")
        logging.info("Use help(function_name) for further documentation")
        logging.info("Use help(module_name) for detailed content list")
        logging.info("")
