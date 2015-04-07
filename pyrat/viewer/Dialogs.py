from PyQt4 import QtGui, QtCore
import pyrat

class FlexInputDialog(QtGui.QDialog):
    def __init__(self, params, title='', doc='no help available', parent=None):
        super(FlexInputDialog, self).__init__(parent)
        self.setWindowTitle(title)
        layout = QtGui.QVBoxLayout(self)

        self.para = params
        self.read = []
        self.doc = doc

        for para in params:
            line = QtGui.QHBoxLayout()
            line.addWidget(QtGui.QLabel(para['text']))
            readline = {'var': para['var'], 'type': para['type'], 'widget': []}
            valuebox = QtGui.QVBoxLayout()
            line.addLayout(valuebox)
            layout.addLayout(line)
            vals = para['value'] if isinstance(para['value'], list) else [para['value']]
            if 'subtext' not in para:
                para.update({'subtext': [' '] * len(vals)})
            hbox = QtGui.QGridLayout()
            valuebox.addLayout(hbox)
            for j, (val, txt) in enumerate(zip(vals, para['subtext'])):
                hbox.addWidget(QtGui.QLabel(txt), j, 0)
                if para['type'] == 'int':
                    wid = QtGui.QSpinBox()
                    if 'range' not in para:
                        para['range'] = [-99999, 99999]
                    wid.setMinimum(para['range'][0])
                    wid.setMaximum(para['range'][1])
                    wid.setValue(val)
                elif para['type'] == 'float':
                    wid = QtGui.QDoubleSpinBox()
                    if 'range' in para:
                        wid.setMinimum(para['range'][0])
                        wid.setMaximum(para['range'][1])
                    wid.setValue(val)
                elif para['type'] == 'bool':
                    wid = QtGui.QCheckBox()
                    if val:
                        wid.setCheckState(2)
                elif para['type'] == 'list':
                    wid = QtGui.QComboBox()
                    wid.addItems(para['range'])
                    wid.setCurrentIndex(para['range'].index(val))
                elif para['type'] in ['openfile', 'opendir', 'savefile']:
                    wid = FlexFilesel(para['type'])
                elif para['type'] == 'str':
                    wid = QtGui.QLineEdit()
                    wid.setText(val)
                readline['widget'].append(wid)
                hbox.addWidget(wid, j, 1)
            # layout.addWidget(HLine())
            self.read.append(readline)

        self.buttons = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel |
                                              QtGui.QDialogButtonBox.Help, QtCore.Qt.Horizontal, self)
        layout.addWidget(self.buttons)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.buttons.helpRequested.connect(self.help)

    def accept(self):
        names = [par['var'] for par in self.para]
        for readline in self.read:
            val = []
            for i, wid in enumerate(readline['widget']):
                if readline['type'] == 'int' or readline['type'] == 'float':
                    val.append(wid.value())
                elif readline['type'] == 'bool':
                    val.append(wid.isChecked())
                elif readline['type'] == 'list':
                    val.append(str(wid.currentText()))
                elif readline['type'] == 'str':
                    val.append(str(wid.text()))
                elif readline['type'] in ['openfile', 'opendir', 'savefile']:
                    val.append(str(wid.text()))
            if len(val) == 1:
                val = val[0]
            # self.para[names.index(readline['name'])]['value'] = val
            # self.para[readline['var']]['value'] = val
            self.para[names.index(readline['var'])]['value'] = val

            setattr(self, readline['var'], val)
        super(FlexInputDialog, self).accept()

    def help(self):
        if self.doc is None:
            self.doc = 'Sorry, no help available'
        try:
            from docutils.core import publish_parts
            doc = publish_parts(self.doc, writer_name='html')['html_body']
        except ImportError:
            doc = self.doc
        foo = QtGui.QMessageBox.information(self, "Routine documentation", doc, QtGui.QMessageBox.Ok)


class FlexFilesel(QtGui.QWidget):
    def __init__(self, type='filename', parent=None):
        self.type = type
        self.value = ''
        super(FlexFilesel, self).__init__(parent)
        layout = QtGui.QHBoxLayout(self)
        self.wid = QtGui.QLineEdit()
        self.wid.setFixedWidth(300)
        layout.addWidget(self.wid)
        self.button = QtGui.QPushButton("Select")
        layout.addWidget(self.button)
        self.button.clicked.connect(self.filesel)

    def filesel(self):
        if self.type == 'openfile':
            self.value = str(QtGui.QFileDialog(self).getOpenFileName())
        elif self.type == 'opendir':
            self.value = str(QtGui.QFileDialog(self).getExistingDirectory())
        elif self.type == 'savefile':
            self.value = str(QtGui.QFileDialog(self).getSaveFileName())
        self.wid.setText(self.value)

    def text(self):
        return self.wid.text()


class PaletteSelector(QtGui.QDialog):
    def __init__(self, tables, current=0, parent=None):
        super(PaletteSelector, self).__init__(parent)
        self.pnames, self.pcolors = tables
        self.palette = current
        self.setWindowTitle("Colour palette chooser")
        layout = QtGui.QVBoxLayout(self)
        self.showpalette = QtGui.QLabel()
        self.update(self.palette)
        wid = QtGui.QLabel("Select colour palette:")
        layout.addWidget(self.showpalette)
        layout.addWidget(wid)
        self.list = QtGui.QListWidget()
        self.list.addItems(self.pnames)
        layout.addWidget(self.list)
        self.buttons = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel, QtCore.Qt.Horizontal, self)
        layout.addWidget(self.buttons)

        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.list.itemClicked.connect(self.item_click)

    def item_click(self, item):
        no = self.pnames.index(item.text())
        self.update(no)

    def update(self, n):
        self.palette = n
        img = QtGui.QImage(self.pcolors[n].tostring(), 256, 1, QtGui.QImage.Format_RGB888)
        self.showpalette.setPixmap(QtGui.QPixmap.fromImage(img).scaled(300,50))


class LayerWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        super(LayerWidget, self).__init__(parent)

        layout = QtGui.QVBoxLayout(self)
        self.treewidget = LayerTreeWidget(parent=self, viewer=parent)
        layout.addWidget(self.treewidget)
        foo = QtGui.QHBoxLayout()
        self.button_bw = QtGui.QRadioButton("B/W mode      ")
        self.button_co = QtGui.QRadioButton("RGB mode      ")
        self.button_bw.setChecked(True)
        foo.addWidget(self.button_bw)
        foo.addWidget(self.button_co)
        layout.addLayout(foo)
        self.button_bw.clicked.connect(self.bwmode)
        self.button_co.clicked.connect(self.rgbmode)
        self.viewer = parent

    def bwmode(self):
        if self.viewer.config['colour']:
            if self.viewer.config['bwlayer'] is not None:
                self.viewer.display[self.viewer.current]['colour'] = False
                self.viewer.processPicture()
            self.redraw()

    def rgbmode(self):
        if self.viewer.config['colour'] is False:
            if not any([foo is None for foo in self.viewer.config['rgblayer']]):
                self.viewer.display[self.viewer.current]['colour'] = True
                self.viewer.processPicture()
            self.redraw()

    def redraw(self):
        if self.viewer.config['colour'] is False:
            self.button_bw.setChecked(True)
        else:
            self.button_co.setChecked(True)
        self.treewidget.redraw()


class LayerTreeWidget(QtGui.QTreeWidget):
    def __init__(self, parent=None, viewer=None):
        QtGui.QTreeWidget.__init__(self, parent)
        self.treelements = {}

        self.viewer = viewer
        self.colors = {
            'default': QtGui.QColor.fromRgb(0, 0, 0, 0),
            'red':   QtGui.QColor.fromRgb(255, 0, 0, 150),
            'green': QtGui.QColor.fromRgb(0, 255, 0, 150),
            'blue':  QtGui.QColor.fromRgb(0, 0, 255, 150),
            'redblue':  QtGui.QColor.fromRgb(255, 0, 255, 150),
            'redgreen':  QtGui.QColor.fromRgb(255, 255, 0, 150),
            'bluegreen':  QtGui.QColor.fromRgb(0, 255, 255, 150),
            'gray':  QtGui.QColor.fromRgb(200, 200, 200, 250)}

        self.setHeaderLabels(['Available Layers'])
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.treeContextMenu)
        self.itemDoubleClicked.connect(self.onDoubleClick)
        self.scalemenu = [None] * 6

    def redraw(self):

        layers = pyrat.data.getLayerNames()
        layers = sorted(layers, key=lambda foo: int(foo.lstrip('/L')))         # sort layers (dict!)
        self.clear()

        for layer in layers:
            lname = layer
            meta = pyrat.data.getAnnotation(layer)

            if 'info' in meta:
                ltext = meta['info']
            else:
                ltext = layer.strip('/')
            self.treelements[lname] = QtGui.QTreeWidgetItem()
            self.treelements[lname].setText(0, ltext)
            self.treelements[lname].setWhatsThis(0, layer)
            self.addTopLevelItem(self.treelements[lname])

            sensor = meta['sensor'] if 'sensor' in meta else 'unknown'
            band = meta['band'] if 'band' in meta else 'unknown'
            query = pyrat.data.queryLayer(layer)
            font = QtGui.QFont("Monospace")
            font.setStyleHint(QtGui.QFont.TypeWriter)
            meta = QtGui.QTreeWidgetItem()
            meta.setText(0, 'meta')
            foo = QtGui.QTreeWidgetItem()
            foo.setFont(0, font)
            foo.setText(0, 'sensor: ' + sensor)
            meta.addChild(foo)
            foo = QtGui.QTreeWidgetItem()
            foo.setFont(0, font)
            foo.setText(0, 'band  : ' + band)
            meta.addChild(foo)
            foo = QtGui.QTreeWidgetItem()
            foo.setFont(0, font)
            foo.setText(0, 'dsize : ' + str(query['dshape']))
            meta.addChild(foo)
            foo = QtGui.QTreeWidgetItem()
            foo.setFont(0, font)
            foo.setText(0, 'lsize : ' + str(query['lshape']))
            meta.addChild(foo)
            foo = QtGui.QTreeWidgetItem()
            foo.setFont(0, font)
            foo.setText(0, 'type  : ' + str(query['dtype']))
            meta.addChild(foo)
            meta.setWhatsThis(0, 'meta')
            self.treelements[lname].addChild(meta)

            channels = pyrat.data.getDataLayerNames(layer)
            metadata = pyrat.data.getAnnotation(layer)
            for k, channel in enumerate(channels):
                cname = channel.split('/')[-1]
                if 'CH_pol' in metadata:
                    cname = metadata['CH_pol'][k]
                self.treelements[channel] = QtGui.QTreeWidgetItem()
                self.treelements[channel].setText(0, cname)
                self.treelements[channel].setWhatsThis(0, channel)
                self.treelements[lname].addChild(self.treelements[channel])
            self.treelements[lname].setExpanded(True)

        self.setFonts()
        self.setColors()

    def setFonts(self):

        activelayer = pyrat.data.active
        if not isinstance(activelayer, list):
            activelayer = [activelayer]

        parent = self.invisibleRootItem()
        for k in range(parent.childCount()):
            litem = parent.child(k)
            font = litem.font(0)
            if litem.whatsThis(0) in activelayer:
                font.setBold(True)
            else:
                font.setBold(False)
            litem.setFont(0, font)

            for l in range(litem.childCount()):
                citem = litem.child(l)
                font = citem.font(0)
                if citem.whatsThis(0) in activelayer:
                    font.setBold(True)
                else:
                    font.setBold(False)
                citem.setFont(0, font)

    def setColors(self):

        parent = self.invisibleRootItem()
        for k in range(parent.childCount()):
            litem = parent.child(k)
            for l in range(litem.childCount()):
                citem = litem.child(l)

                if self.viewer.config['colour'] is False:
                    if citem.whatsThis(0) == self.viewer.config['bwlayer']:
                        citem.setBackgroundColor(0, self.colors['gray'])
                    else:
                        citem.setBackground(0, self.colors['default'])
                else:
                    if citem.whatsThis(0) == self.viewer.config['rgblayer'][0] and \
                            citem.whatsThis(0) == self.viewer.config['rgblayer'][1] and \
                            citem.whatsThis(0) == self.viewer.config['rgblayer'][2]:
                        citem.setBackgroundColor(0, self.colors['gray'])
                    elif citem.whatsThis(0) == self.viewer.config['rgblayer'][0] and \
                            citem.whatsThis(0) == self.viewer.config['rgblayer'][1]:
                        citem.setBackgroundColor(0, self.colors['redgreen'])
                    elif citem.whatsThis(0) == self.viewer.config['rgblayer'][0] and \
                            citem.whatsThis(0) == self.viewer.config['rgblayer'][2]:
                        citem.setBackgroundColor(0, self.colors['redblue'])
                    elif citem.whatsThis(0) == self.viewer.config['rgblayer'][1] and \
                            citem.whatsThis(0) == self.viewer.config['rgblayer'][2]:
                        citem.setBackgroundColor(0, self.colors['bluegreen'])
                    elif citem.whatsThis(0) == self.viewer.config['rgblayer'][0]:
                        citem.setBackgroundColor(0, self.colors['red'])
                    elif citem.whatsThis(0) == self.viewer.config['rgblayer'][1]:
                        citem.setBackgroundColor(0, self.colors['green'])
                    elif citem.whatsThis(0) == self.viewer.config['rgblayer'][2]:
                        citem.setBackgroundColor(0, self.colors['blue'])
                    else:
                        citem.setBackground(0, self.colors['default'])

    def treeContextMenu(self, position):
        item = self.currentItem()
        itemname = item.whatsThis(0).split('/')[-1]
        menu = QtGui.QMenu()
        if itemname[0] is 'L':
            itemname = item.whatsThis(0)
            menu.addAction("show / activate", lambda: self.activate(itemname))
            menu.addSeparator()
            menu.addAction("add to activation", lambda: self.addactive(itemname))
            menu.addAction("remove activation", lambda: self.delactive(itemname))
            menu.addSeparator()
            menu.addAction("delete data set", lambda: self.delayer(itemname))
            menu.addSeparator()
            scalemenu = menu.addMenu("scaling method")

            scaling = self.viewer.display[itemname]['scaling']

            methods = ['amplitude', 'intensity', 'phase', '0.0->1.0', 'min->max', 'lables']
            for meth in methods:
                foo = QtGui.QAction(meth, scalemenu, checkable=True)
                if scaling == meth:
                    foo.setChecked(True)
                foo.changed.connect(lambda m=meth: self.scaling(itemname, m))
                scalemenu.addAction(foo)

        elif itemname[0] is 'D':
            itemname = item.whatsThis(0)
            itemlayer = '/' + itemname.split('/')[1]

            menu.addAction("add to activation", lambda: self.addactive(itemname))
            menu.addAction("remove activation", lambda: self.delactive(itemname))
            menu.addSeparator()
            if self.viewer.config['colour'] is True:
                viewlayer = '/' + self.viewer.config['rgblayer'][0].split('/')[1]
            else:
                viewlayer = '/' + self.viewer.config['bwlayer'].split('/')[1]
            if itemlayer == viewlayer:
                if self.viewer.config['colour'] is True:
                    menu.addAction("set layer as red ", lambda: self.setRGB(itemname, 'r'))
                    menu.addAction("set layer as green ", lambda: self.setRGB(itemname, 'g'))
                    menu.addAction("set layer as blue ", lambda: self.setRGB(itemname, 'b'))
                else:
                    menu.addAction("show layer", lambda: self.setRGB(itemname))

        menu.exec_(self.viewport().mapToGlobal(position))

    def setRGB(self, channel, col=None):
        if col == 'r':
            self.viewer.config['rgblayer'][0] = channel
        elif col == 'g':
            self.viewer.config['rgblayer'][1] = channel
        elif col == 'b':
            self.viewer.config['rgblayer'][2] = channel
        else:
            self.viewer.config['bwlayer'] = channel
        self.redraw()
        self.viewer.processPicture()

    def delayer(self, layer):
        pyrat.data.delLayer(layer)
        self.viewer.updateViewer()
        self.redraw()
        self.setFonts()

    def activate(self, layer):
        print("activate:", layer)
        pyrat.data.activateLayer(layer)
        self.redraw()
        if 'D' in layer:
            lay = '/'+layer.split('/')[1]
            self.viewer.display[lay]['bwlayer'] = layer
            self.viewer.display[lay]['colour'] = False
        self.viewer.updateViewer(layer=layer)

    def addactive(self, layer):
        activelayer = pyrat.data.active
        if not isinstance(activelayer, list):
            activelayer = [activelayer]
        if layer not in activelayer:
            activelayer.append(layer)
            pyrat.data.activateLayer(activelayer)
            self.redraw()
            self.setFonts()

    def delactive(self, layer):
        activelayer = pyrat.data.active
        if not isinstance(activelayer, list):
            activelayer = [activelayer]
        if layer in activelayer and len(activelayer) > 1:
            activelayer.remove(layer)
            pyrat.data.activateLayer(activelayer)
            self.redraw()
            self.setFonts()
        if len(activelayer) == 1:
            self.viewer.updateViewer(layer=activelayer[0])

    def scaling(self, layer, method):
        oldmethod = self.viewer.display[layer]['scaling']
        if method != oldmethod:
            self.viewer.display[layer]['scaling'] = method
            if self.viewer.config['colour'] is True:
                viewlayer = '/' + self.viewer.config['rgblayer'][0].split('/')[1]
            else:
                viewlayer = '/' + self.viewer.config['bwlayer'].split('/')[1]
            if viewlayer == layer:
                self.viewer.showCurrentLayer(force=True)
            else:
                GenPyramid(layer=layer, force=True).run()

    def onDoubleClick(self, item):
        layer = item.whatsThis(0)
        print("ondoubleclick", layer)
        if layer[1] == 'L':
            self.activate(layer)
