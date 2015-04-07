from PyQt4 import QtGui, QtCore


class CropBoxWidget(QtGui.QWidget):
    """
    Custom widget for crop region queries (4 values)
    """
    def __init__(self, title=None, parent=None):
        super(CropBoxWidget, self).__init__(parent)
        self.value = [0]*4
        layout = QtGui.QGridLayout(self)
        if isinstance(title, str):
            layout.addWidget(QtGui.QLabel(title), 0, 0)
        self.crop1 = QtGui.QSpinBox()
        self.crop2 = QtGui.QSpinBox()
        self.crop3 = QtGui.QSpinBox()
        self.crop4 = QtGui.QSpinBox()
        cropbox = QtGui.QGridLayout()
        cropbox.addWidget(QtGui.QLabel("range start"), 0, 0)
        cropbox.addWidget(self.crop1, 1, 0)
        cropbox.addWidget(QtGui.QLabel("range end"), 2, 0)
        cropbox.addWidget(self.crop2, 3, 0)
        cropbox.addWidget(QtGui.QLabel("azimuth start"), 0, 1)
        cropbox.addWidget(self.crop3, 1, 1)
        cropbox.addWidget(QtGui.QLabel("azimuth end"), 2, 1)
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


class FileselWidget(QtGui.QWidget):
    """
    Custom widget for file selection
    """
    def __init__(self, type='openfile', title=None, parent=None):
        super(FileselWidget, self).__init__(parent)
        self.type = type
        self.value = ''
        mainlayout = QtGui.QGridLayout(self)
        if isinstance(title, str):
            mainlayout.addWidget(QtGui.QLabel(title), 0, 0)
        layout = QtGui.QHBoxLayout()
        self.text = QtGui.QLineEdit()
        self.text.setFixedWidth(300)
        layout.addWidget(self.text)
        self.button = QtGui.QPushButton("Select")
        layout.addWidget(self.button)
        mainlayout.addLayout(layout, 1, 0)
        self.button.clicked.connect(self.filesel)
        self.setContentsMargins(0, 0, 0, 0)
        self.layout().setContentsMargins(0, 0, 0, 0)

    def filesel(self):
        if self.type == 'openfile':
            self.value = str(QtGui.QFileDialog(self).getOpenFileName())
        elif self.type == 'opendir':
            self.value = str(QtGui.QFileDialog(self).getExistingDirectory())
        elif self.type == 'savefile':
            self.value = str(QtGui.QFileDialog(self).getSaveFileName())

        self.text.setText(self.value)

    def setvalue(self, value):
        self.value = value
        self.text.setText(value)

    def getvalue(self):
        return self.text.text()


class ProductContentWidget(QtGui.QWidget):
    """
    Custom widget for product component selection
    """
    def __init__(self, title=None, parent=None, bands=True, polar=True, products=None):
        super(ProductContentWidget, self).__init__(parent)
        self.bandflag = bands
        self.polflag = polar
        mainlayout = QtGui.QGridLayout(self)
        if isinstance(title, str):
            mainlayout.addWidget(QtGui.QLabel(title), 0, 0)
        layout = QtGui.QGridLayout()
        layout.addWidget(QtGui.QLabel("Product"), 0, 0)
        self.product = QtGui.QComboBox()
        self.product.addItems(products)
        layout.addWidget(self.product, 1, 0)
        if self.bandflag is True:
            layout.addWidget(QtGui.QLabel("Band"), 0, 1)
            self.band = QtGui.QComboBox()
            self.band.addItems(["*"])
            layout.addWidget(self.band, 1, 1)
        if self.polflag is True:
            layout.addWidget(QtGui.QLabel("Polarisation"), 0, 2)
            self.polar = QtGui.QComboBox()
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


class HLine(QtGui.QFrame):
    """
    Class for drawining a horiyontal line in a qt widget
    """
    def __init__(self, parent=None):
        super(HLine, self).__init__(parent)
        self.setFrameShape(QtGui.QFrame.HLine)
        self.setFrameShadow(QtGui.QFrame.Sunken)


class VLine(QtGui.QFrame):
    """
    Class for drawining a vertical line in a qt widget
    """
    def __init__(self, parent=None):
        super(VLine, self).__init__(parent)
        self.setFrameShape(QtGui.QFrame.VLine)
        self.setFrameShadow(QtGui.QFrame.Sunken)


class DelayedUpdater(QtGui.QWidget):

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
