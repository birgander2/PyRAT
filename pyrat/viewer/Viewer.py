import pyrat
import numpy as np
import logging
from PyQt5 import QtGui, QtCore, QtWidgets
from .Dialogs import PaletteSelector, LayerWidget
from .StatusBar import *
from . import egg
from pyrat.tools import multimap

from pyrat.tools import colortables
from pyrat.viewer.tools import sarscale, subsample


class MainWindow(QtWidgets.QMainWindow):

    modules = []

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setWindowTitle("PyRAT - Radar Tools")

        self.box = [0, 100, 0, 100]  # current display image coordinates
        self.rubberbox = [0, 1, 0, 1]  # coordinates of a rubberbox selection

        self.size = [100, 100]
        self.factor = 1.0
        self.sarscale = 2.5
        self.type = 'A'

        self.data = []  # link to preview data
        self.current = None  # currently displayed layer
        self.display = {}  # dictionary of display configs
        self.config = {}  # current config

        self.block_redraw = False
        self.picture_redraw = False
        self.show_rubberband = False
        self.undolist = []

        self.makeActions()
        self.makeToolbar()
        self.makeStatusBar()
        self.makeMenu()
        self.makeView()
        self.initPlugins()
        self.resize(1000, 800)
        self.central.setSizes([150, 800])
        self.updateDisplayList()
        self.dragX = 0
        self.dragY = 0
        self.imgwidth = 0
        self.imgheight = 0
        self.show()
        self.rubberband = QtWidgets.QRubberBand(QtWidgets.QRubberBand.Rectangle, self.imageLabel)
        self.palette = 0

        QtWidgets.QShortcut(QtGui.QKeySequence("Ctrl+T"), self, self.easterEgg)
        self.central.setMinimumSize(100, 100)

    def makeToolbar(self):
        self.openTB = QtWidgets.QAction(QtGui.QIcon('icons/document-open.png'), 'Open', self,
                                        triggered=lambda: pyrat.load.RatHDF.guirun(self))
        # self.closeTB = QtWidgets.QAction(QtGui.QIcon('icons/document-close.png'), 'Close', self)
        # self.zoominTB = QtGui.QAction(QtGui.QIcon('icons/zoom-in.png'), 'Zoom in', self)
        # self.zoomoutTB = QtGui.QAction(QtGui.QIcon('icons/zoom-out.png'), 'Zoom out', self)
        # self.zoomresetTB = QtGui.QAction(QtGui.QIcon('icons/zoom-fit-best.png'), 'Fit zoom', self)
        self.seperatorTB = QtWidgets.QAction(self)

        self.toolbar1 = self.addToolBar("File")
        self.toolbar1.addAction(self.openTB)
        # self.toolbar1.addAction(self.closeTB)

        self.toolbar2 = self.addToolBar("Display")
        self.toolbar2.addAction(self.zoomOutAct)
        self.viewCombo = QtWidgets.QComboBox(self)
        self.viewCombo.insertItems(1, ["100%", "Fit to window", "Fit to width", "Fit to height", "100%"])
        self.viewCombo.setEditable(False)
        self.viewCombo.activated.connect(self.comboZoom)
        self.toolbar2.addWidget(self.viewCombo)
        self.toolbar2.addAction(self.zoomInAct)
        self.toolbar2.addAction(self.paletteAct)

        # self.toolbar3 = self.addToolBar("Layer")
        # self.toolbar3.addAction(self.zoominTB)
        # self.toolbar3.addAction(self.zoomresetTB)

    # -------------------------------- STATUS BAR
    def makeStatusBar(self):
        self.statusBar = StatusBar(self)
        self.statusBar.setMessage(message='no data loaded')

    # -------------------------------- DEFAULT ACTIONS (only those not implemented trough plugins)
    def makeActions(self):
        self.exitAct = QtWidgets.QAction('Exit', self, shortcut='Q', triggered=self.close)
        self.zoomInAct = QtWidgets.QAction(QtGui.QIcon('icons/zoom-in.png'), "Zoom &In (25%)", self, shortcut="up",
                                           triggered=lambda: self.zoom(3.0 / 2.0))
        self.zoomOutAct = QtWidgets.QAction(QtGui.QIcon('icons/zoom-out.png'), "Zoom &Out (25%)", self, shortcut="down",
                                            triggered=lambda: self.zoom(2.0 / 3.0))
        self.fitToWindowAct = QtWidgets.QAction(QtGui.QIcon('icons/zoom-fit-best.png'), "Reset view", self,
                                                shortcut="f",
                                                triggered=self.resetView)
        self.viewAmpAct = QtWidgets.QAction("View as amplitude", self, checkable=True, shortcut="1",
                                            triggered=self.viewAsAmplitude)
        self.viewPhaAct = QtWidgets.QAction("View as phase", self, checkable=True, shortcut="2",
                                            triggered=self.viewAsPhase)
        self.viewCohAct = QtWidgets.QAction("View as coherence", self, checkable=True, shortcut="3",
                                            triggered=self.viewAsCoherence)
        self.viewBrighter = QtWidgets.QAction("View brighter", self, shortcut="right", triggered=self.brighterView)
        self.viewDarker = QtWidgets.QAction("View darker", self, shortcut="left", triggered=self.darkerView)
        self.undoAct = QtWidgets.QAction('Undo', self, shortcut='Ctrl+z', triggered=self.undo)
        self.paletteAct = QtWidgets.QAction(QtGui.QIcon('icons/color_wheel.png'), 'Palette', self,
                                            triggered=self.paletteChooser)

    def easterEgg(self):
        tetris = egg.Tetris(self)
        tetris.show()

    # ---------------------------------- PLUGINS
    def initPlugins(self):
        from inspect import getmembers, isclass

        modules = self.modules.copy()
        modules += [pyrat.filter, pyrat.load, pyrat.save, pyrat.transform,
                    pyrat.polar, pyrat.insar, pyrat.plugins, pyrat.viewer]

        logging.debug("Scanning for GUI elements:")
        for current_module in modules:
            modules = getmembers(current_module, isclass)
            for mod in modules:
                if issubclass(mod[1], pyrat.Worker):
                    # plugin = mod[1]()
                    plugin = mod[1]
                    if hasattr(plugin, 'gui'):
                        logging.debug(" Attaching GUI element : " + mod[0])
                        plugin.registerGUI(self)

    def makeMenu(self):
        self.menubar = self.menuBar()
        self.menue = {
            "File": self.menubar.addMenu('File'),
            "General": self.menubar.addMenu('General'),
            "Tools": self.menubar.addMenu('Tools'),
            "SAR": self.menubar.addMenu('SAR'),
            "PolSAR": self.menubar.addMenu('PolSAR'),
            "InSAR": self.menubar.addMenu('InSAR')}

        if 'oss' not in pyrat.__version__:
            self.menue["DLR"] = self.menubar.addMenu('DLR')

        self.menue["Help"] = self.menubar.addMenu('Help')

        self.menue.update({"File|Import raster": self.menue["File"].addMenu('Import raster')})
        self.menue.update({"File|Import spaceborne": self.menue["File"].addMenu('Import spaceborne')})
        self.menue.update({"File|Import airborne": self.menue["File"].addMenu('Import airborne')})
        self.menue.update({"File|Import pixmap": self.menue["File"].addMenu('Import pixmap')})
        foo = self.menue["File"].addSeparator()
        foo.setWhatsThis("File|line1")
        self.menue.update({"File|Export to raster": self.menue["File"].addMenu('Export to raster')})
        self.menue.update({"File|Export to pixmap": self.menue["File"].addMenu('Export to pixmap')})
        foo = self.menue["File"].addSeparator()
        foo.setWhatsThis("File|line2")
        self.menue["File"].addAction(self.exitAct)

        self.menue["General"].addAction(self.undoAct)
        self.menue["General"].addSeparator()
        self.menue["General"].addAction(self.fitToWindowAct)
        self.menue["General"].addAction(self.zoomInAct)
        self.menue["General"].addAction(self.zoomOutAct)
        self.menue["General"].addSeparator()
        self.menue["General"].addAction(self.viewBrighter)
        self.menue["General"].addAction(self.viewDarker)
        self.menue["General"].addSeparator()
        self.viewSel = QtWidgets.QActionGroup(self.menue["General"])
        foo = self.viewSel.addAction(self.viewAmpAct)
        self.menue["General"].addAction(foo)
        foo = self.viewSel.addAction(self.viewPhaAct)
        self.menue["General"].addAction(foo)
        foo = self.viewSel.addAction(self.viewCohAct)
        self.menue["General"].addAction(foo)

        self.menue.update({"SAR|Speckle filter": self.menue["SAR"].addMenu('Speckle filter')})
        self.menue.update({"SAR|Edge detection": self.menue["SAR"].addMenu('Edge detection')})
        self.menue.update({"SAR|Texture": self.menue["SAR"].addMenu('Texture')})
        self.menue.update({"SAR|line1": self.menue["SAR"].addSeparator()})
        self.menue.update({"SAR|Geometry": self.menue["SAR"].addMenu('Geometry')})
        self.menue.update({"SAR|Amplitude": self.menue["SAR"].addMenu('Amplitude')})
        self.menue.update({"SAR|line2": self.menue["SAR"].addSeparator()})
        self.menue.update({"SAR|Sidelobe control": self.menue["SAR"].addMenu('Sidelobe control')})
        self.menue.update({"SAR|Spectral tools": self.menue["SAR"].addMenu('Spectral tools')})
        self.menue.update({"SAR|line3": self.menue["SAR"].addSeparator()})
        self.menue.update({"SAR|Transform": self.menue["SAR"].addMenu('Transform')})
        self.menue.update({"PolSAR|Calibration": self.menue["PolSAR"].addMenu('Calibration')})
        self.menue.update({"PolSAR|Speckle filter": self.menue["PolSAR"].addMenu('Speckle filter')})
        self.menue.update({"PolSAR|Change detection": self.menue["PolSAR"].addMenu('Change detection')})
        self.menue.update({"PolSAR|Classification": self.menue["PolSAR"].addMenu('Classification')})
        self.menue.update({"PolSAR|Decompositions": self.menue["PolSAR"].addMenu('Decompositions')})
        self.menue.update({"PolSAR|Parameters": self.menue["PolSAR"].addMenu('Parameters')})
        self.menue.update({"PolSAR|Transform": self.menue["PolSAR"].addMenu('Transform')})
        self.menue.update({"InSAR|Phase noise filter": self.menue["InSAR"].addMenu('Phase noise filter')})
        self.menue.update({"InSAR|Transform": self.menue["InSAR"].addMenu('Transform')})
        self.menue.update({"Tools|Geometry": self.menue["Tools"].addMenu('Geometry')})
        self.menue.update({"Tools|Filter": self.menue["Tools"].addMenu('Filter')})

        # self.menue["View"].addAction(self.viewAmpAct)
        # self.menue["View"].addAction(self.viewPhaAct)
        # self.menue["View"].addAction(self.viewCohAct)

    # ---------------------------------- VIEW AREA
    def makeView(self):
        # self.central = QtGui.QWidget()
        # self.HLayout = QtGui.QHBoxLayout(self.central)
        self.central = QtWidgets.QSplitter(self)
        self.central.setOpaqueResize(False)
        self.central.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        self.tree = LayerWidget(self)
        # self.tree.setFixedWidth(200)
        # self.tree.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        self.central.addWidget(self.tree)

        # self.HLayout.addItem(self.spacer)

        self.frame = QtWidgets.QWidget()
        # self.frame = DelayedUpdater()

        self.frame.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.imageLabel = QtWidgets.QLabel(self.frame)
        self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)
        self.imageLabel.setStyleSheet("QLabel { background-color: #333 }")
        self.imageLabel.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.imageLabel.setScaledContents(False)
        self.imageLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.imageLabel.resize(self.frame.width(), self.frame.height())
        self.central.addWidget(self.frame)

        self.setCentralWidget(self.central)

        self.resizeframeevent = self.frame.resizeEvent
        self.frame.resizeEvent = self.resizeFrameEvent  # a bit dirtyt to overload like this...

        # self.central.move.connect(lambda: self.splitterMoved())

        # self.frame = QtGui.QWidget()
        # self.imageLabel = QtGui.QLabel(self.frame)
        # self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)
        # self.imageLabel.setStyleSheet("QLabel { background-color: #333 }")
        # self.imageLabel.setSizePolicy(QtGui.QSizePolicy.Ignored, QtGui.QSizePolicy.Ignored)
        # self.imageLabel.setScaledContents(False)
        # self.imageLabel.resize(self.frame.width(), self.frame.height())
        # self.setCentralWidget(self.frame)

    # def split(self, event):
    #     # QtGui.QWidget.resizeEvent(self.frame, event)
    #     # self.central.update()
    #     self.resizeEvent(event)
    #     # self.central.update()
    #     print("Event", self.frame.width(), self.frame.height())

    def paletteChooser(self):
        wid = PaletteSelector(colortables(), current=0, parent=self)
        res = wid.exec_()
        if res == 1:
            self.display[self.current]['palette'] = wid.palette
            self.processPicture()

    def updateDisplayList(self):
        layers = pyrat.data.getLayerIDs()
        if not isinstance(layers, list):
            layers = [layers]

        for layer in layers:
            if layer not in self.display:
                config = {'colour': None, 'bwlayer': None, 'rgblayer': [None, None, None],
                          'palette': 0, 'scaling': 'amplitude'}
                dsets = pyrat.data.getDataLayerIDs(layer)
                if len(dsets) == 1:
                    config['bwlayer'] = dsets[0]
                    config['colour'] = False
                elif len(dsets) == 2:
                    config['bwlayer'] = dsets[0]
                    config['rgblayer'] = [dsets[0], dsets[1], dsets[1]]
                    config['colour'] = True
                else:
                    config['bwlayer'] = dsets[0]
                    config['rgblayer'] = [dsets[0], dsets[2], dsets[1]]
                    config['colour'] = True
                self.display.update({layer: config})

    def updateViewer(self, method=None, layer=None):

        if layer is None:
            current = pyrat.data.active
        else:
            current = layer
        if isinstance(current, list):
            current = current[-1]
        if current is None:
            if len(pyrat.data.layers) != 0:  # show first layer if existing
                current = list(pyrat.data.layers.keys())[0]
            else:
                if hasattr(self, 'data'):
                    del self.data
                self.current = None
                self.imageLabel.clear()
        if current is not None:
            self.undolist.append((list(pyrat.data.getLayerIDs()), pyrat.data.active))
            self.current = '/' + current.split('/')[1]
            self.updateDisplayList()
            self.config = self.display[self.current]
            if method is not None:
                self.config['scaling'] = method
            self.tree.redraw()
            self.showCurrentLayer()

    def showCurrentLayer(self, force=False):
        layer = self.current
        self.statusBar.setMessage(message=' Updating view ', colour='R')
        self.getPyramid(layer, force=force)
        self.processPicture(fitwin=1)
        self.statusBar.setMessage(message=' Ready ', colour='G')

    def getPyramid(self, layer, force=False):
        mode = self.display[layer]['scaling']
        self.data = GenPyramid(layer=layer, force=force, mode=mode).run()
        self.size = self.data[0].shape[-2:][::-1]  # [xmin,xmax,ymin,ymax]

    # -----------------------------------------------------------
    # -----------------------------------------------------------
    # -----------------------------------------------------------
    # -----------------------------------------------------------

    def data2img(self, cut_box, scale=0):
        if self.config['colour'] is False:
            channel = [int(self.config['bwlayer'].split('D')[1])]
        else:
            channel = [int(foo.split('D')[1]) for foo in self.config['rgblayer']]
        method = self.config['scaling']

        img = self.data[scale][..., cut_box[2]:cut_box[3], cut_box[0]:cut_box[1]][channel, ...].copy()

        if img.dtype == 'uint8':
            img = img.copy()
        elif method == 'amplitude' or method == 'intensity' or method == None:
            if method == 'intensity':
                img = np.abs(img) ** 0.35
            else:
                img = np.abs(img) ** 0.7
            if self.config['colour'] is True:
                for k in range(img.shape[0]):
                    img[k, ...] /= np.mean(img[k, ...])
            img = sarscale(img, factor=self.sarscale)
        elif method == 'phase':
            img = np.uint8(np.clip(img / np.pi * 128 + 127, 0.0, 255.0))
        elif method == '0.0->1.0':
            img = np.uint8(np.clip(img * 255.0, 0.0, 255.0))
        elif method == 'min->max':
            img -= np.min(img)
            img = np.uint8(img / np.max(img) * 255.0)
        elif method == 'labels':
            img -= np.min(self.data[0])
            img = np.uint8(img / np.max(self.data[0]) * 255.0)
        else:
            img = np.uint8(int)

        # img = img[..., 0:img.shape[-2] // 4 * 4, 0:img.shape[-1] // 4 * 4]    # QT Limitation!!
        img = np.rollaxis(np.rollaxis(img, axis=2), axis=2)

        if self.config['palette'] != 0:
            self.img = colortables(self.config['palette'])[1][img]
        else:
            self.img = img

        # img = np.rot90(np.rollaxis(np.rollaxis(img, axis=2), axis=2))

        if self.config['colour'] is True:
            return QtGui.QImage(img.tostring(), img.shape[1], img.shape[0], QtGui.QImage.Format_RGB888)
        else:
            return QtGui.QImage(img.tostring(), img.shape[1], img.shape[0], QtGui.QImage.Format_Indexed8)

    def processPicture(self, **kwargs):
        if "fitwin" in kwargs.keys():  # self.size  = size of data set
            self.scale = len(self.data) - 1  # self.scale = pyramid level
            while self.scale > 0:  # self.box   = displayed part
                if self.data[self.scale].shape[-1] > self.imageLabel.width() or \
                                self.data[self.scale].shape[-2] > self.imageLabel.height():
                    break
                self.scale -= 1
            self.box = [0, self.size[0], 0, self.size[1]]
            cut_box = [foo // 2 ** self.scale for foo in self.box]
        else:
            self.scale = len(self.data) - 1
            while self.scale > 0:
                cut_box = [foo // 2 ** self.scale for foo in self.box]
                if cut_box[1] - cut_box[0] > self.imageLabel.width() or \
                                        cut_box[3] - cut_box[2] > self.imageLabel.height():
                    break
                self.scale -= 1

            cut_box = [foo // 2 ** self.scale for foo in self.box]

        cut_box[1] = cut_box[0] + (cut_box[1] - cut_box[0]) // 4 * 4  # QT Limitation!!
        cut_box[3] = cut_box[2] + (cut_box[3] - cut_box[2]) // 4 * 4
        self.scale = int(np.clip(self.scale, 0, len(self.data) - 1))
        self.box[1] = cut_box[1] * int(2 ** self.scale)
        self.box[3] = cut_box[3] * int(2 ** self.scale)

        img = self.data2img(cut_box, scale=self.scale)

        xWin = self.imageLabel.width()
        yWin = self.imageLabel.height()
        winRatio = 1.0 * xWin / yWin

        self.imgwidth = img.width()
        self.imgheight = img.height()

        imgRatio = 1.0 * self.imgwidth / self.imgheight

        if imgRatio >= winRatio:  # match widths
            self.imgwidth = xWin
            self.imgheight = int(xWin / imgRatio)
        else:  # match heights
            self.imgheight = yWin
            self.imgwidth = int(yWin * imgRatio)

        self.factor = int(100.0 * self.imgwidth / (self.box[1] - self.box[0]))
        if self.factor <= 100:
            img = img.scaled(self.imgwidth, self.imgheight)  # Bilinear?
        else:
            img = img.scaled(self.imgwidth, self.imgheight)  # Nearest Neighbour

        self.statusBar.setMessage(size=1, zoom=1, level=1, scale=1)
        self.viewCombo.setItemText(0, str(int(self.factor)) + '%')

        # colortable = [QtGui.qRgb(i, i, i) for i in range(256)]
        p = colortables(self.config['palette'])[1]
        colortable = [QtGui.qRgb(p[i, 0], p[i, 1], p[i, 2]) for i in range(256)]

        img.setColorTable(colortable)
        self.imageLabel.setPixmap(QtGui.QPixmap.fromImage(img))

    def zoom(self, factor, mx=0, my=0):
        px = self.box[0] + int(
            (1.0 * self.box[1] - self.box[0]) / self.imageLabel.width() * (mx + self.imageLabel.width() // 2))
        py = self.box[2] + int(
            (1.0 * self.box[3] - self.box[2]) / self.imageLabel.height() * (my + self.imageLabel.height() // 2))

        midx = self.box[0] + (self.box[1] - self.box[0]) // 2
        midy = self.box[2] + (self.box[3] - self.box[2]) // 2
        sizx = self.box[1] - self.box[0]
        sizy = self.box[3] - self.box[2]
        newx = np.clip(int(sizx / factor), 0, self.size[0])
        newy = np.clip(int(sizy / factor), 0, self.size[1])
        imgRatio = 1.0 * newx / newy
        xWin = self.imageLabel.width()
        yWin = self.imageLabel.height()
        winRatio = 1.0 * xWin / yWin
        if imgRatio >= winRatio:  # match widths:
            newy = np.clip(int(newx / winRatio), 0, self.size[1])
        else:
            newx = np.clip(int(newy * winRatio), 0, self.size[0])
        newx = np.clip(newx, 8, self.size[0])
        newy = np.clip(newy, 8, self.size[1])
        if midx - newx // 2 < 0:
            midx = newx // 2
        if midx + newx // 2 > self.size[0]:
            midx = self.size[0] - newx // 2
        if midy - newy // 2 < 0:
            midy = newy // 2
        if midy + newy // 2 > self.size[1]:
            midy = self.size[1] - newy // 2
        self.box = [midx - newx // 2, midx + newx // 2, midy - newy // 2, midy + newy // 2]

        if mx != 0 and my != 0:
            midx = px - int((1.0 * self.box[1] - self.box[0]) / self.imageLabel.width() * mx)
            midy = py - int((1.0 * self.box[3] - self.box[2]) / self.imageLabel.height() * my)
            sizx = self.box[1] - self.box[0]
            sizy = self.box[3] - self.box[2]
            if midx - sizx // 2 < 0:
                midx = sizx // 2
            if midx + sizx // 2 > self.size[0]:
                midx = self.size[0] - sizx // 2
            if midy - sizy // 2 < 0:
                midy = sizy // 2
            if midy + sizy // 2 > self.size[1]:
                midy = self.size[1] - sizy // 2
            self.box = [midx - sizx // 2, midx + sizx // 2, midy - sizy // 2, midy + sizy // 2]

        if self.picture_redraw is False:  # needed to avoid concurrent redraws
            self.picture_redraw = True
            if hasattr(self, 'data'):
                self.processPicture()
            self.picture_redraw = False

    def undo(self):
        if len(self.undolist) >= 2:
            actual = set(self.undolist[-1][0])
            before = set(self.undolist[-2][0])
            diff = actual.difference(before)
            pyrat.data.activateLayer(self.undolist[-2][1])
            pyrat.data.delLayer(list(diff))
            for layer in list(diff):
                del self.display[layer]
            self.undolist.pop()
            self.undolist.pop()
            self.updateViewer()

    def viewAsAmplitude(self):
        self.type = 'A'
        if hasattr(self, 'data'):
            self.processPicture()

    def viewAsCoherence(self):
        self.type = 'C'
        if hasattr(self, 'data'):
            self.processPicture()

    def viewAsPhase(self):
        self.type = 'P'
        if hasattr(self, 'data'):
            self.processPicture()

    def darkerView(self):
        self.sarscale += 0.5
        if hasattr(self, 'data'):
            self.processPicture()

    def brighterView(self):
        self.sarscale -= 0.5
        if self.sarscale < 0.5:
            self.sarscale = 0.5
        if hasattr(self, 'data'):
            self.processPicture()

    def resetView(self):
        self.sarscale = 2.5
        if hasattr(self, 'data'):
            self.processPicture(fitwin=1)

    def comboZoom(self, index):
        if index == 1:
            if hasattr(self, 'data'):
                self.processPicture(fitwin=1)
            self.viewCombo.setCurrentIndex(0)
        elif index == 2:
            print("Not implemented")
            self.viewCombo.setCurrentIndex(0)
        elif index == 3:
            print("Not implemented")
            self.viewCombo.setCurrentIndex(0)
        elif index == 4:
            print("Not implemented")
            self.viewCombo.setCurrentIndex(0)

    def getPosVal(self, x, y):
        """
        Calc position of mouse in data coordinates and returns string containing the
        values at that position.
        """
        hscale = (self.box[1] - self.box[0]) / self.frame.rect().width()
        vscale = (self.box[3] - self.box[2]) / self.frame.rect().height()
        scale = max(hscale, vscale)

        dpx = (self.imageLabel.width() - self.imgwidth) // 2
        dpy = (self.imageLabel.height() - self.imgheight) // 2

        posx = self.box[0] + int((x - dpx) * scale)
        posy = self.box[2] + int((y - dpy) * scale)
        if self.current is None:
            return None
        if 0 <= posx < self.size[0] and 0 <= posy < self.size[1]:
            values = pyrat.data.getData(block=(posy, posy + 1, posx, posx + 1), layer=self.current)
            if values.shape == ():
                values = values.reshape(1)
            txt = '<pre>Cursor position: [' + str(posy) + ', ' + str(posx) + ']'
            for k, val in enumerate(values):
                txt += '<br>D' + str(k) + ':   '
                if np.iscomplexobj(val):
                    txt += str(np.abs(val)) + ' abs  /  ' + str(np.angle(val, deg=True)) + ' deg'
                elif self.config['scaling'] == 'phase':
                    txt += str(val * 180.0 / np.pi) + ' deg'
                else:
                    txt += str(val)
            txt += '</pre>'
        else:
            txt = ''
        return txt

    def resizeFrameEvent(self, event):
        if self.block_redraw is False:  # needed to avoid concurrent redraws -> segfault
            self.block_redraw = True
            self.imageLabel.setGeometry(0, 0, self.frame.width(), self.frame.height())
            midx = self.box[0] + (self.box[1] - self.box[0]) // 2
            midy = self.box[2] + (self.box[3] - self.box[2]) // 2
            sizx = self.box[1] - self.box[0]
            sizy = self.box[3] - self.box[2]
            newx = np.clip(int(sizx), 0, self.size[0])
            newy = np.clip(int(sizy), 0, self.size[1])
            imgRatio = 1.0 * newx / newy
            xWin = self.imageLabel.width()
            yWin = self.imageLabel.height()

            winRatio = 1.0 * xWin / yWin
            if imgRatio >= winRatio:  # match widths:
                newy = np.clip(int(newx / winRatio), 0, self.size[1])
            else:
                newx = np.clip(int(newy * winRatio), 0, self.size[0])
            midx -= (midx - newx // 2) * ((midx - newx // 2) < 0)
            midy -= (midy - newy // 2) * ((midy - newy // 2) < 0)

            self.box = [midx - newx // 2, midx + newx // 2, midy - newy // 2, midy + newy // 2]

            if hasattr(self, 'data') and len(self.data) > 0:
                self.processPicture()
            self.resizeframeevent(event)
            self.block_redraw = False

    def wheelEvent(self, event):
        x = event.x() - self.central.x() - self.frame.x()  # easier with self.frame.mapFrom()
        y = event.y() - self.central.y() - self.frame.y()
        if hasattr(self, 'data') and self.current is not None:
            if event.angleDelta().y() < 0:
                self.zoom(2.0 / 3.0, mx=x - self.imageLabel.width() // 2,
                          my=y - self.imageLabel.height() // 2)
            if event.angleDelta().y() > 0:
                self.zoom(3.0 / 2.0, mx=x - self.imageLabel.width() // 2,
                          my=y - self.imageLabel.height() // 2)
            self.processPicture()

    def mousePressEvent(self, event):
        if self.show_rubberband is True:
            # self.origin = event.pos()
            self.origin = self.frame.mapFrom(self, event.pos())
            self.rubberband.setGeometry(QtCore.QRect(self.origin, QtCore.QSize()))
            self.rubberband.show()
        else:
            if event.button() == QtCore.Qt.RightButton and self.current is not None:
                xoff = self.central.x() + self.frame.x()  # easier with self.frame.mapFrom()
                yoff = self.central.y() + self.frame.y()
                x = event.x() - xoff
                y = event.y() - yoff
                txt = self.getPosVal(x, y)
                if hasattr(self, "posval"):
                    self.posval.setText(txt)
                    self.posval.adjustSize()
                    self.posval.show()
                else:
                    self.posval = QtWidgets.QLabel(txt)
                    self.posval.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
                    self.posval.setWindowTitle("position / value")
                    self.posval.setGeometry(self.x() + xoff, self.y() + yoff, 200, 100)
                    self.posval.adjustSize()
                    self.posval.show()
            else:
                self.dragX = event.x()
                self.dragY = event.y()

        QtWidgets.QWidget.mousePressEvent(self, event)

    def mouseReleaseEvent(self, event):
        if self.show_rubberband is True:
            if self.rubberband.isVisible():
                self.rubberband.hide()
                self.show_rubberband = False

                crop = self.rubberband.geometry()

                wx = self.imageLabel.width()
                wy = self.imageLabel.height()
                win_ratio = wx / wy

                ix = self.box[1] - self.box[0]
                iy = self.box[3] - self.box[2]
                im_ratio = ix / iy

                if im_ratio >= win_ratio:  # width matches
                    scale = ix / wx
                else:
                    scale = iy / wy

                xs = int(ix / scale)
                ys = int(iy / scale)
                xo = (wx - xs) // 2
                yo = (wy - ys) // 2

                x1 = self.box[0] + int(scale * (crop.x() - xo))
                x2 = x1 + int(scale * crop.width())
                y1 = self.box[2] + int(scale * (crop.y() - yo))
                y2 = y1 + int(scale * crop.height())

                x1 = np.clip(x1, self.box[0], self.box[1])
                x2 = np.clip(x2, self.box[0], self.box[1])
                y1 = np.clip(y1, self.box[2], self.box[3])
                y2 = np.clip(y2, self.box[2], self.box[3])

                self.rubberbox = [y1, y2, x1, x2]

        else:
            if event.button() == QtCore.Qt.LeftButton:
                dx = int((1.0 * self.box[1] - self.box[0]) / self.imageLabel.width() * (self.dragX - event.x()))
                dy = int((1.0 * self.box[3] - self.box[2]) / self.imageLabel.height() * (self.dragY - event.y()))
                midx = self.box[0] + (self.box[1] - self.box[0]) // 2 + dx
                midy = self.box[2] + (self.box[3] - self.box[2]) // 2 + dy
                sizx = self.box[1] - self.box[0]
                sizy = self.box[3] - self.box[2]
                if midx - sizx // 2 < 0:
                    midx = sizx // 2
                if midx + sizx // 2 > self.size[0]:
                    midx = self.size[0] - sizx // 2
                if midy - sizy // 2 < 0:
                    midy = sizy // 2
                if midy + sizy // 2 > self.size[1]:
                    midy = self.size[1] - sizy // 2
                self.box = [midx - sizx // 2, midx + sizx // 2, midy - sizy // 2, midy + sizy // 2]
                if hasattr(self, 'data') and len(self.data) > 0:
                    self.processPicture()
                self.imageLabel.setGeometry(0, 0, self.frame.width(), self.frame.height())
            elif event.button() == QtCore.Qt.RightButton:
                # if hasattr(self, "posval"):
                #     self.posval.hide()
                pass
        QtWidgets.QWidget.mouseReleaseEvent(self, event)

    def mouseMoveEvent(self, event):
        if self.show_rubberband is True:
            if self.rubberband.isVisible():
                pos = self.frame.mapFrom(self, event.pos())
                # todo: limit size
                self.rubberband.setGeometry(QtCore.QRect(self.origin, pos).normalized())
        else:
            if event.buttons() == QtCore.Qt.LeftButton:
                self.imageLabel.move(event.x() - self.dragX, event.y() - self.dragY)
            elif event.buttons() == QtCore.Qt.RightButton:
                if hasattr(self, "posval"):
                    xoff = self.central.x() + self.frame.x()
                    yoff = self.central.y() + self.frame.y()
                    x = event.x() - xoff
                    y = event.y() - yoff
                    txt = self.getPosVal(x, y)
                    self.posval.setText(txt)
                    self.posval.adjustSize()
        QtWidgets.QWidget.mouseMoveEvent(self, event)


class GenPyramid():
    """
    Generation of the presumming pyramid

    :author: Andreas Reigber
    """

    def __init__(self, layer=None, force=False, mode='amplitude', *args, **kwargs):
        super(GenPyramid, self).__init__(*args, **kwargs)
        if layer is None:
            self.layer = pyrat.data.active
        else:
            self.layer = layer
        self.hdfgroup = pyrat.data.layers[self.layer].group

        query = pyrat.data.queryLayer(self.layer)
        self.force = force  # force recalculation
        self.mode = mode  # rescaling method
        self.scale = 0
        self.lshape = (np.prod(query['lshape']),)
        self.dshape = query['dshape']
        self.cut = [dim // 2 * 2 for dim in self.dshape]
        self.shp = [dim // 2 * 2 for dim in self.dshape]
        self.dset = []

        self.nblock = 0
        self.progress = pyrat.tools.ProgressBar('  Updating View', 3 * self.dshape[-2] / 128)

    def __del__(self):
        del self.progress

    def run(self, *args, **kwargs):

        # print("pyramid level", self.scale)
        if self.scale == 0:
            ilay = 'D'
            olay = 'P/0'
            ishp = [dim // 2 * 2 for dim in self.dshape]
            oshp = [dim // 2 * 2 for dim in self.dshape]
        else:
            ilay = 'P/' + str(self.scale - 1)
            olay = 'P/' + str(self.scale)
            ishp = [dim // 2 * 2 for dim in self.dshape]
            oshp = [dim // 2 for dim in self.dshape]

        idat = self.hdfgroup[ilay]
        if olay in self.hdfgroup and self.force is False:
            self.dset.append(self.hdfgroup[olay])
            self.scale += 1
            self.dshape = oshp
            self.progress.update(self.scale * 100)
            self.run()
        else:
            odat = self.hdfgroup.require_dataset(olay, self.lshape + tuple(oshp), 'float32')
            self.dset.append(odat)
            if min(oshp) > 1:
                idx, ivalid, ibs = self.calc_blocks(ishp[-2], pack=True, blocksize=128)
                odx, ovalid, obs = self.calc_blocks(oshp[-2], pack=False, blocksize=ibs // (ishp[-2] // oshp[-2]))
                nb = 0
                for bidx in idx:
                    inputs = []

                    for ix in bidx:
                        data = idat[..., ix[0]:ix[1], 0:ishp[-1]]
                        inputs.append((subsample, (data, self.lshape + (obs, oshp[1]), self.mode)))
                    result = multimap(inputs, mode='method')
                    # result = map(absrebin, inputs)
                    for res in result:
                        odat[..., odx[nb][0] + ovalid[nb][0]:odx[nb][1], :] = res[..., ovalid[nb][0]:ovalid[nb][1], :]
                        nb += 1
                        self.nblock += 1
                        self.progress.update(self.nblock)
                self.scale += 1
                self.dshape = oshp
                self.run()
        return self.dset

    def calc_blocks(self, size, blocksize=128, pack=False):

        nthreads = pyrat._nthreads
        if blocksize > size:  # maximum equal image size
            blocksize = size

        blocks = [[0, blocksize]]
        while blocks[-1][1] < size:
            blocks.append([blocks[-1][1], blocks[-1][1] + blocksize])
        offset = blocks[-1][1] - size  # last block starts earlier
        blocks[-1][0] -= offset  # with increased overlap
        blocks[-1][1] -= offset

        valid = [0] * len(blocks)  # calculate the valid part of each block (start, end)
        for k, block in enumerate(blocks):
            if k == len(blocks) - 1:  # last block
                if k == 0:
                    valid[k] = block  # only one block
                else:
                    valid[k] = [blocks[-2][1] - block[0], block[1] - block[0]]
            else:  # middle block
                valid[k] = [0, block[1] - block[0]]

        if pack is True:
            if len(blocks) > 1 and nthreads > 1:
                blocks = [blocks[i:i + nthreads] for i in range(0, len(blocks), nthreads)]
            else:
                blocks = [[block] for block in blocks]
        return blocks, valid, blocksize
