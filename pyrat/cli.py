from __future__ import print_function


def info():
    from pyrat import data

    data.info()


def activate(layer, silent=False):
    from pyrat import data

    data.activateLayer(layer, silent=silent)


def delete(layer, silent=False):
    from pyrat import data

    data.delLayer(layer, silent=silent)


def getdata(**kwargs):
    from pyrat import data

    return data.getData(**kwargs)


def adddata(array):
    from pyrat import data

    return data.addLayer(array)


def getmeta(**kwargs):
    from pyrat import data

    return data.getAnnotation(**kwargs)


def setmeta(meta, **kwargs):
    from pyrat import data

    data.setAnnotation(meta, **kwargs)


def show():
    import pyrat
    import sys
    from PyQt4 import QtGui

# TODO: Currently only the last active data set can be viewed!
    if len(pyrat.data.active) > 0:
        active = pyrat.data.active[0]
        arr = pyrat.data.getData()
        if isinstance(arr, list):
            arr = arr[0]


        if arr.ndim == 3:
            arr = arr[0:3, ...]
        if arr.ndim == 4:
            arr = arr[0:3, ...]
    app = QtGui.QApplication(sys.argv)
    pyrat.app = pyrat.viewer.MainWindow()
    pyrat.app.updateViewer()
    sys.exit(app.exec_())

gui = show

# def plist():
#     import sys
#     from inspect import getmembers, isclass
#
#     current_module = sys.modules[__name__]
#     modules = getmembers(current_module, isclass)
#     print()
#     print("Content of module " + __name__ + ":")
#     print()
#     for mod in modules:
#         if 'pyrat' in mod[1].__module__:
#             doc = str(mod[1].__doc__)
#             if doc != 'None':
#                 doc = doc.split('\n')[1]
#             print(mod[0].ljust(20) + doc)
