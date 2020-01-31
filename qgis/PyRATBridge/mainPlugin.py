"""
PyRATBridge
===========

This Plugin imports the functionality of PyRAT into QGIS

:author: Felix Weinmann <felix.weinmann@dlr.de>

"""

from qgis.core import QgsTask, QgsTaskManager, Qgis, QgsProject
from qgis.PyQt.QtWidgets import QAction, QFileDialog, QInputDialog, QDockWidget
from qgis.PyQt.QtCore import Qt
from qgis.utils import iface
import copy
import numpy as np
from os import path
try:
    import pyrat
    from pyrat.viewer.Dialogs import FlexInputDialog, LayerTreeWidget
    pyratImport = True
except ImportError:
    pyratImport = False

qgis_types = {
    0: None,  # UnknownDataType
    1: "int8",
    2: "uint16",
    3: "int16",
    4: "uint32",
    5: "int32",
    6: "float32",
    7: "float64",
    8: None,  # CInt16
    9: None,  # CInt32
    10: "complex32",
    11: "complex64",
    12: None,  # ARGB32. Color, alpha, red, green, blue
    13: None,  # ARGB32_Premultiplied alpha, red, green, blue
}


class PyRATBridge:
    """This is the main plugin class for GUI and the connection to PyRAT"""

    def __init__(self):
        self.taskManager = QgsTaskManager()
        if pyratImport:
            pyrat.viewer.GenPyramid = GenPyramidInterface
            pyrat.tools.ProgressBar = ProgressBarToQGIS

    def unload(self):
        """Cleanup when disabling the plugin"""
        if pyratImport:
            PyRATBridge.clearPyRAT()
            self.pyratMenu.clear()
            iface.removeDockWidget(self.pyratLayerTree)
            self.taskManager.cancelAll()
            ViewerToQGISInterface.display.clear()

    def addMenuEntry(self, pyratTool):
        """Adds a PyRAT Tool to the QGIS-Menu"""
        menus = pyratTool.gui['menu'].split('|')
        submenu = self.pyratMenu
        for menu in menus:
            if menu not in [action.text() for action in submenu.actions()]:
                submenu = submenu.addMenu(menu)
            else:
                submenu = [action.menu() for action in submenu.actions() if
                           action.text() == menu][0]

        action = QAction(pyratTool.gui['entry'], iface.mainWindow())
        action.triggered.connect(lambda:
                                 PyRATBridge.menuAction(self, pyratTool))
        submenu.addAction(action)

    def initGui(self):
        """Initalise the Plugin-UI"""
        if not pyratImport:
            iface.messageBar().pushMessage("PyRAT not found!",
                                           level=Qgis.Critical)
            return

        if 'PyRAT' not in [action.text() for action in
                           iface.mainWindow().menuBar().actions()]:
            self.pyratMenu = iface.mainWindow().menuBar().addMenu('PyRAT')
        else:
            self.pyratMenu = [action.menu() for action in
                              iface.mainWindow().menuBar().actions() if
                              action.text() == 'PyRAT'][0]

        action = QAction("Layer2PyRAT", iface.mainWindow())
        action.triggered.connect(PyRATBridge.layerToPyrat)
        self.pyratMenu.addAction(action)

        action = QAction("PyRAT2Layer", iface.mainWindow())
        action.triggered.connect(PyRATBridge.pyratToLayer)
        self.pyratMenu.addAction(action)

        action = QAction("Cleanup PyRAT", iface.mainWindow())
        action.triggered.connect(PyRATBridge.clearPyRAT)
        self.pyratMenu.addAction(action)

        action = QAction("Show PyRAT GUI", iface.mainWindow())
        action.triggered.connect(self.showPyrat)
        self.pyratMenu.addAction(action)

        self.pyratMenu.addSeparator()

        # Init PyRAT-Tools, adapted from pyrat.viewer for qgis
        from inspect import getmembers, isclass

        modules = [pyrat.load, pyrat.save, pyrat.transform, pyrat.filter,
                   pyrat.polar, pyrat.insar, pyrat.plugins, pyrat.viewer]

        for current_module in modules:
            modules = getmembers(current_module, isclass)
            for mod in modules:
                if issubclass(mod[1], pyrat.Worker):
                    plugin = mod[1]
                    if(hasattr(plugin, 'gui') and
                       plugin.gui['entry'] != "Python console"):
                        self.addMenuEntry(plugin)

        self.pyratLayerTree = QDockWidget("PyRAT Layers", iface.mainWindow())
        PyRATBridge.layerTreeWidget = LayerTreeWidget(
                parent=self.pyratLayerTree,
                viewer=ViewerToQGISInterface)
        self.pyratLayerTree.setObjectName("PyRAT Layers")
        self.pyratLayerTree.setWidget(PyRATBridge.layerTreeWidget)
        iface.addDockWidget(Qt.LeftDockWidgetArea, self.pyratLayerTree)

    def menuAction(self, pyratTool):
        """Start pyratTool after Menu-Click"""
        para_backup = copy.deepcopy(pyratTool.para)

        if 'name' not in dir(pyratTool):
            pyratTool.name = pyratTool.__name__

        if len(pyratTool.para) > 0:
            dlg = FlexInputDialog(pyratTool.para,
                                  parent=iface.mainWindow(),
                                  title=pyratTool.name,
                                  doc=pyratTool.__doc__)

        if len(pyratTool.para) == 0 or dlg.exec_() == 1:
            task = PyRATTask(pyratTool, para_backup)
            self.taskManager.addTask(task)

    def layerToPyrat():
        """Imports a QGIS-Layer into PyRAT"""
        layers = list()
        for layer in QgsProject.instance().layerTreeRoot().layerOrder():
            # 1: QgsMapLayer.LayerType.RasterLayer
            if layer.type() == 1:
                layers.append(layer.name())

        layername, s = QInputDialog.getItem(
                iface.mainWindow(),
                "Select a layer",
                "Select a layer to export to PyRAT:",
                layers,
                editable=False)
        if not s:
            return

        layer = QgsProject.instance().mapLayersByName(layername)[0]
        dataProv = layer.dataProvider()
        extent = dataProv.extent()
        rows = layer.height()
        cols = layer.width()
        block = dataProv.block(1, extent, cols, rows)
        arr = np.frombuffer(block.data(),
                            dtype=qgis_types[block.dataType()]
                            ).reshape((rows, cols))
        pyratlayer = pyrat.adddata(arr)

        # Add metadata to the PyRAT-Layer
        description = layer.crs().description()
        meta = {"info": layer.name(),
                "geo_min_east": extent.xMinimum(),
                # Subtract 1 due to QGIS inclusive minimum
                "geo_min_north": extent.yMinimum() - 1,
                "geo_ps_east": layer.rasterUnitsPerPixelX(),
                "geo_ps_north": layer.rasterUnitsPerPixelY()}

        if description.startswith('WGS 84 / UTM zone '):
            zone = int(description[:-1].rsplit(" ", 1)[1])
            if description[-1] == "S":
                zone = -zone
            meta["geo_projection"] = 1
            meta["geo_zone"] = zone

        pyrat.setmeta(meta)
        ViewerToQGISInterface.display[pyratlayer] = {'scaling': 'min->max',
                                                     'bwlayer': pyratlayer,
                                                     'colour': False}
        PyRATBridge.layerTreeWidget.redraw()

    def pyratToLayer(layerid=None):
        """Exports a PyRAT-layer into QGIS"""
        if type(layerid) is str:
            pyrat.data.activateLayer(layerid)
        annotation = pyrat.data.getAnnotation()
        if 'info' in annotation:
            filename = path.join(pyrat.data.tmpdir, annotation['info'] +
                                 ".rat")
        else:
            filename = path.join(pyrat.data.tmpdir, "PyRAT.rat")

        filename, s = QFileDialog.getSaveFileName(
                iface.mainWindow(),
                "Save the PyRAT-Layer",
                filename,
                "RAT-File (*.rat)")

        if not s or filename == "":
            return

        pyrat.save.rat((filename, "rat"), geo_envi_hdr=True)
        iface.addRasterLayer(filename, path.basename(filename).split(".")[0])

    def showPyrat(self):
        pyrat.show()

    def clearPyRAT():
        pyrat.pyrat_reset()
        ViewerToQGISInterface.display.clear()
        PyRATBridge.layerTreeWidget.redraw()


class ViewerToQGISInterface:
    """This Class is a 'viewer' for pyrats LayerTree Widget shown in QGIS"""

    config = {'colour': False, 'bwlayer': "/Undefined",
              'rgblayer': (None, None, None)}
    display = {}

    def updateViewer(layer=None):
        pass


class GenPyramidInterface:
    """
    This class replaces pyrat.viewer.GenPyramid to disable
    the scaling method options in the LayerTree Widget in QGIS
    """

    def __init__(self, layer=None, force=None, mode=None):
        pass

    def run(self):
        pass


class ProgressBarToQGIS:
    """
    Disables the ProgressBar to prevent crashes with opened QGIS Python Console
    """

    def __init__(self, message, max, width=None):
        pass

    def __del__(self):
        pass

    def update(self, val):
        pass


class PyRATTask(QgsTask):
    """This class handles the async execution of a PyRAT-Tool"""

    def __init__(self, pyratTool, para_backup):
        QgsTask.__init__(self)
        self.pyratTool = pyratTool
        self.para_backup = para_backup
        self.failed = False
        self.guionly = False
        self.layer = None
        self.existinglayers = list()

    def run(self):
        """The async executed code"""
        self.plugin = self.pyratTool()
        self.plugin.crash_handler = self.crashHandler
        self.existinglayers = pyrat.data.getLayerIDs()
        self.layer = self.plugin.run()

        setattr(self.pyratTool, 'para', self.para_backup)
        return self.layer is not False

    def crashHandler(self, ex):
        """
        Overrides the PyRAT crash handler to prevent
        the termination of QGIS
        """
        try:
            raise ex
        except AttributeError:
            # Gui-only Plugins
            self.guionly = True
        except Exception:
            self.failed = True
            raise ex

    def finished(self, result):
        """
        This function is threadsafe for GUI-Actions and
        called after run terminates.
        """
        if self.guionly:
            self.pyratTool.guirun(iface.mainWindow())

        if result and not self.failed:
            iface.messageBar().pushMessage(self.pyratTool.name + " finished.",
                                           level=Qgis.Success)

            for layer in [newlayer for newlayer in pyrat.data.getLayerIDs()
                          if newlayer not in self.existinglayers]:
                # Show the generated Layer(s) in QGIS
                anno = pyrat.data.getAnnotation(layer=layer)
                if 'info' not in anno:
                    anno['info'] = "Pyrat-Layer " + layer
                pyrat.data.setAnnotation({'info': anno['info'] + "-" +
                                          self.pyratTool.name},
                                         layer=layer)
                ViewerToQGISInterface.display[layer] = {'scaling': 'min->max',
                                                        'bwlayer': layer,
                                                        'colour': False}
                PyRATBridge.pyratToLayer(self.layer)
            PyRATBridge.layerTreeWidget.redraw()
        else:
            iface.messageBar().pushMessage(self.pyratTool.name +
                                           " failed. Look in the (system)" +
                                           " console for more information.",
                                           level=Qgis.Critical)
        del self.plugin
