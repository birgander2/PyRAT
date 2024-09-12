from PyQt5 import QtCore, QtWidgets


class CropBoxWidget(QtWidgets.QWidget):
    """
    Custom widget for crop region queries (4 values)
    """
    def __init__(self, title=None, parent=None):
        super(CropBoxWidget, self).__init__(parent)
        self.value = [0]*4
        layout = QtWidgets.QGridLayout(self)
        if isinstance(title, str):
            layout.addWidget(QtWidgets.QLabel(title), 0, 0)
        self.crop1 = QtWidgets.QSpinBox()
        self.crop2 = QtWidgets.QSpinBox()
        self.crop3 = QtWidgets.QSpinBox()
        self.crop4 = QtWidgets.QSpinBox()
        cropbox = QtWidgets.QGridLayout()
        cropbox.addWidget(QtWidgets.QLabel("range start"), 0, 0)
        cropbox.addWidget(self.crop1, 1, 0)
        cropbox.addWidget(QtWidgets.QLabel("range end"), 2, 0)
        cropbox.addWidget(self.crop2, 3, 0)
        cropbox.addWidget(QtWidgets.QLabel("azimuth start"), 0, 1)
        cropbox.addWidget(self.crop3, 1, 1)
        cropbox.addWidget(QtWidgets.QLabel("azimuth end"), 2, 1)
        cropbox.addWidget(self.crop4, 3, 1)
        layout.addLayout(cropbox, 1, 0)
        self.minmaxMemory = []
        self.setrange([[0, 9999999], [0, 9999999], [0, 9999999], [0, 9999999]])
        self.setContentsMargins(0, 0, 0, 0)
        self.layout().setContentsMargins(0, 0, 0, 0)

    def setvalues(self, value):
        self.crop1.setValue(value[2])
        self.crop2.setValue(value[3])
        self.crop3.setValue(value[0])
        self.crop4.setValue(value[1])

    def getvalues(self):
        self.value[2] = self.crop1.value()
        self.value[3] = self.crop2.value()
        self.value[0] = self.crop3.value()
        self.value[1] = self.crop4.value()
        return self.value

    def setrange(self, minmax):
        #function is called two times
        #print(minmax)
        #[[0, 9999999], [0, 9999999], [0, 9999999], [0, 9999999]]
        #[[0, 0], [0, 0], [0, 0], [0, 0]]
        # saves the para from the first call in self.minmaxMemory
        if len(self.minmaxMemory)<= 0:
            self.minmaxMemory = minmax
        else:
            minmax = self.minmaxMemory

        self.crop1.setMinimum(minmax[2][0])
        self.crop1.setMaximum(minmax[2][1])
        self.crop2.setMinimum(minmax[3][0])
        self.crop2.setMaximum(minmax[3][1])
        self.crop3.setMinimum(minmax[0][0])
        self.crop3.setMaximum(minmax[0][1])
        self.crop4.setMinimum(minmax[1][0])
        self.crop4.setMaximum(minmax[1][1])


class FileselWidget(QtWidgets.QWidget):
    """
    Custom widget for file selection
    """
    def __init__(self, type='openfile', title=None, parent=None):
        super(FileselWidget, self).__init__(parent)
        self.type = type
        self.value = ''
        mainlayout = QtWidgets.QGridLayout(self)
        if isinstance(title, str):
            mainlayout.addWidget(QtWidgets.QLabel(title), 0, 0)
        layout = QtWidgets.QHBoxLayout()
        self.text = QtWidgets.QLineEdit()
        self.text.setFixedWidth(300)
        layout.addWidget(self.text)
        self.button = QtWidgets.QPushButton("Select")
        layout.addWidget(self.button)
        mainlayout.addLayout(layout, 1, 0)
        self.button.clicked.connect(self.filesel)
        self.setContentsMargins(0, 0, 0, 0)
        self.layout().setContentsMargins(0, 0, 0, 0)

    def filesel(self):
        if self.type == 'openfile':
            # self.value = str(QtWidgets.QFileDialog(self).getOpenFileName()[0])
            self.value = str(QtWidgets.QFileDialog.getOpenFileName(parent=self,
                                                                   options=QtWidgets.QFileDialog.DontUseNativeDialog)[0])
        elif self.type == 'opendir':
            self.value = str(QtWidgets.QFileDialog.getExistingDirectory(parent=self,
                                                                       options=QtWidgets.QFileDialog.DontUseNativeDialog))
            # self.value = str(QtWidgets.QFileDialog(self).getExistingDirectory())
        elif self.type == 'savefile':
            self.value = str(QtWidgets.QFileDialog().getSaveFileName(parent=self,
                                                                     options=QtWidgets.QFileDialog.DontUseNativeDialog)[0])

        self.text.setText(self.value)

    def setvalue(self, value):
        self.value = value
        self.text.setText(value)

    def getvalue(self):
        return self.text.text()


class ProductContentWidget(QtWidgets.QWidget):
    """
    Custom widget for product component selection
    """
    def __init__(self, title=None, parent=None, bands=True, polar=True, products=None):
        super(ProductContentWidget, self).__init__(parent)
        self.bandflag = bands
        self.polflag = polar
        mainlayout = QtWidgets.QGridLayout(self)
        if isinstance(title, str):
            mainlayout.addWidget(QtWidgets.QLabel(title), 0, 0)
        layout = QtWidgets.QGridLayout()
        layout.addWidget(QtWidgets.QLabel("Product"), 0, 0)
        self.product = QtWidgets.QComboBox()
        self.product.addItems(products)
        layout.addWidget(self.product, 1, 0)
        if self.bandflag is True:
            layout.addWidget(QtWidgets.QLabel("Band"), 0, 1)
            self.band = QtWidgets.QComboBox()
            self.band.addItems(["*"])
            layout.addWidget(self.band, 1, 1)
        if self.polflag is True:
            layout.addWidget(QtWidgets.QLabel("Polarisation"), 0, 2)
            self.polar = QtWidgets.QComboBox()
            self.polar.addItems(["*"])
            layout.addWidget(self.polar, 1, 2)
        mainlayout.addLayout(layout, 1, 0)
        self.setContentsMargins(0, 0, 0, 0)
        self.layout().setContentsMargins(0, 0, 0, 0)

    def getvalue(self, index):
        if index == 0:
            return str(self.product.currentText())
        elif index == 1 and self.bandflag is True:
            return str(self.band.currentText())
        elif index == 2 and self.polflag is True:
            return str(self.polar.currentText())

    def setvalue(self, index, val):
        if index == 0:
            self.product.setCurrentIndex(self.product.findText(val))
        elif index == 1 and self.bandflag is True:
            self.band.setCurrentIndex(self.band.findText(val))
        elif index == 2 and self.polflag is True:
            self.polar.setCurrentIndex(self.polar.findText(val))

    def updatepolar(self, pols):
        # for k in range(self.polar.count()):
        #     self.polar.removeItem(k)
        self.polar.clear()
        self.polar.addItem("*")
        for p in pols:
            self.polar.addItem(p)

    def updatebands(self, bands):
        # for k in range(self.band.count()):
        #     self.band.removeItem(k)
        self.band.clear()
        self.band.addItem("*")
        for b in bands:
            self.band.addItem(b)

    def updateproducts(self, products):
        self.product.clear()
        for p in products:
            self.product.addItem(p)


class ProductContentWidget_UAVSAR(QtWidgets.QWidget):
    """
    Custom widget for product component selection UAVSAR
    """
    def __init__(self, title=None, parent=None, bands=True, polar=True, products=None, tracks=None):
        super(ProductContentWidget_UAVSAR, self).__init__(parent)
        self.bandflag = bands
        self.polflag = polar
        mainlayout = QtWidgets.QGridLayout(self)
        if isinstance(title, str):
            mainlayout.addWidget(QtWidgets.QLabel(title), 0, 0)
        layout = QtWidgets.QGridLayout()
        layout.addWidget(QtWidgets.QLabel("Product"), 0, 0)
        self.product = QtWidgets.QComboBox()
        self.product.addItems(products)
        layout.addWidget(self.product, 1, 0)

        layout.addWidget(QtWidgets.QLabel("Tracks"), 0, 1)
        self.tracks = QtWidgets.QComboBox()
        self.tracks.addItems(tracks)
        layout.addWidget(self.tracks, 1, 1)
        if self.bandflag is True:
            layout.addWidget(QtWidgets.QLabel("Band"), 0, 2)
            self.band = QtWidgets.QComboBox()
            self.band.addItems(["*"])
            layout.addWidget(self.band, 1, 2)
        if self.polflag is True:
            layout.addWidget(QtWidgets.QLabel("Polarisation"), 0, 3)
            self.polar = QtWidgets.QComboBox()
            self.polar.addItems(["*"])
            layout.addWidget(self.polar, 1, 3)
        mainlayout.addLayout(layout, 1, 0)
        self.setContentsMargins(0, 0, 0, 0)
        self.layout().setContentsMargins(0, 0, 0, 0)

    def getvalue(self, index):
        if index == 0:
            return str(self.product.currentText())
        elif index == 1:
            return str(self.tracks.currentText())
        elif index == 2 and self.bandflag is True:
            return str(self.band.currentText())
        elif index == 3 and self.polflag is True:
            return str(self.polar.currentText())

    def setvalue(self, index, val):
        if index == 0:
            self.product.setCurrentIndex(self.product.findText(val))
        elif index == 1:
            self.tracks.setCurrentIndex(self.tracks.findText(val))
        elif index == 2 and self.bandflag is True:
            self.band.setCurrentIndex(self.band.findText(val))
        elif index == 3 and self.polflag is True:
            self.polar.setCurrentIndex(self.polar.findText(val))

    def updatepolar(self, pols):
        # for k in range(self.polar.count()):
        #     self.polar.removeItem(k)
        self.polar.clear()
        self.polar.addItem("*")
        for p in pols:
            self.polar.addItem(p)

    def updatebands(self, bands):
        # for k in range(self.band.count()):
        #     self.band.removeItem(k)
        self.band.clear()
        self.band.addItem("*")
        for b in bands:
            self.band.addItem(b)

    def updateproducts(self, products):
        self.product.clear()
        for p in products:
            self.product.addItem(p)

    def updatetracks(self, tracks):
        self.tracks.clear()
        self.tracks.addItem("*")
        for p in tracks:
            self.tracks.addItem(p)


class BoolWidget(QtWidgets.QWidget):
    def __init__(self, text=None, parent=None):
        super(BoolWidget, self).__init__(parent)
        mainlayout = QtWidgets.QHBoxLayout(self)
        self.text = QtWidgets.QLabel(text)
        mainlayout.addWidget(self.text, 0)
        self.checkbox = QtWidgets.QCheckBox()
        mainlayout.addWidget(self.checkbox, 1)

    def getvalue(self):
        return self.checkbox.isChecked()


class HLine(QtWidgets.QFrame):
    """
    Class for drawining a horiyontal line in a qt widget
    """
    def __init__(self, parent=None):
        super(HLine, self).__init__(parent)
        self.setFrameShape(QtWidgets.QFrame.HLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class VLine(QtWidgets.QFrame):
    """
    Class for drawining a vertical line in a qt widget
    """
    def __init__(self, parent=None):
        super(VLine, self).__init__(parent)
        self.setFrameShape(QtWidgets.QFrame.VLine)
        self.setFrameShadow(QtWidgets.QFrame.Sunken)


class DelayedUpdater(QtWidgets.QWidget):

    def __init__(self):
        super(DelayedUpdater, self).__init__()
        self.delayEnabled = True
        self.delayTimeout = 500
        self._resizeTimer = QtCore.QTimer(self)
        self._resizeTimer.timeout.connect(self._delayedUpdate)

    def resizeEvent(self, event):
        if self.delayEnabled:
            self._resizeTimer.start(self.delayTimeout)
            self.setUpdatesEnabled(False)
        super(DelayedUpdater, self).resizeEvent(event)

    def _delayedUpdate(self):
        self._resizeTimer.stop()
        self.setUpdatesEnabled(True)
