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



class ViewerTree(QtGui.QTreeWidget):

    def __init__(self, parent=None):
        QtGui.QTreeWidget.__init__(self, parent)
        self.setHeaderLabels(['Layer     ', 'Content'])
        self.setColumnWidth(0, 75)
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.treeContextMenu)
        self.itemDoubleClicked.connect(self.onDoubleClick)

    def addLayers(self, elements):
        parent = self.invisibleRootItem()
        for layer, channel in elements.items():
            entry = self.addParent(parent, 0, layer.strip('/'), 'what parent?')
            for ch in channel:
                self.addChild(entry, 0, ch.lstrip(layer), 'what child?')

    def addParent(self, parent, column, title, data):
        item = QtGui.QTreeWidgetItem(parent, [title])
        item.setData(column, QtCore.Qt.UserRole, data)
        item.setChildIndicatorPolicy(QtGui.QTreeWidgetItem.ShowIndicator)
        item.setExpanded (True)
        # item.setCheckState (column, QtCore.Qt.Unchecked)
        return item

    def addChild(self, parent, column, title, data):
        item = QtGui.QTreeWidgetItem(parent, [title])
        item.setData(column, QtCore.Qt.UserRole, data)
        #item.setCheckState (column, QtCore.Qt.Unchecked)
        return item

    def delTree(self):
        """
        Deletes all layers in tree
        """
        root = self.invisibleRootItem()
        self.clear()

    def activateRow(self, layers):
        if isinstance(layers, str):
            layers = [layers]
        root = self.invisibleRootItem()
        child_count = root.childCount()
        for i in range(child_count):
            item = root.child(i)
            layer = '/'+str(item.text(0))
            if layer in layers:
                item.setBackground(0, QtGui.QBrush(QtCore.Qt.lightGray))
                item.setBackground(1, QtGui.QBrush(QtCore.Qt.lightGray))
            else:
                item.setBackground(0, QtGui.QBrush(QtCore.Qt.white))
                item.setBackground(1, QtGui.QBrush(QtCore.Qt.white))

    def treeContextMenu(self, position):
        item = self.currentItem()
        if str(item.text(0))[0] is 'L':
            menu = QtGui.QMenu()
            layer = '/'+str(item.text(0))
            menu.addAction("Activate layer", lambda: pyrat.app.activateLayer(layer))
            menu.addAction("Delete layer ", lambda: pyrat.app.delLayer(layer))
            menu.exec_(self.viewport().mapToGlobal(position))
        elif str(item.text(0))[0] is 'D':
            pass

    # def handleChanged(self, item, column):
    #     if item.checkState(column) == QtCore.Qt.Checked:
    #         print "checked", item, item.text(column)
    #     if item.checkState(column) == QtCore.Qt.Unchecked:
    #         print "unchecked", item, item.text(column)
    #

    def onDoubleClick(self, *item):
        """
        param:
        item[0] = QTreeWidgetItem
        item[1] = integer
        """
        if item[0].childCount() > 0:
            layer = '/'+str(item[0].text(0))
            pyrat.app.activateLayer(layer)
        elif item[0].parent().childCount() > 0:
            layer = '/'+str(item[0].parent().text(0))
            pyrat.app.activateLayer(layer)


