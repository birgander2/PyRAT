from PyQt5 import QtCore, QtGui, QtWidgets


class StatusBar:
    def __init__(self, parent):
        self.parent = parent
        self.statusBar = parent.statusBar()


        self.pixmap_green  = QtGui.QPixmap("icons/green.png", "PNG")
        self.pixmap_yellow = QtGui.QPixmap("icons/yellow.png", "PNG")
        self.pixmap_red    = QtGui.QPixmap("icons/red.png", "PNG")
        self.statusPixmap = QtWidgets.QLabel("image", self.statusBar)
        self.statusPixmap.setPixmap(self.pixmap_green)
        
        self.__statusTimer = QtCore.QTimer(self.parent)
        # self.parent.connect(self.__statusTimer, QtCore.SIGNAL("timeout()"), self.resetMessage)
        self.__statusLabel = QtWidgets.QLabel("Default", self.statusBar)
        self.progressbar = QtWidgets.QProgressBar(self.statusBar)
        self.progressbar.setMinimum(0)
        self.progressbar.setMaximum(100)
        self.progressbar.setValue(0)
        self.progressbar.setTextVisible(False)

        self.lastMessage = ''
       
        self.statusSize  = QtWidgets.QLabel("", self.statusBar)
        self.statusZoom  = QtWidgets.QLabel("", self.statusBar)
        self.statusScale = QtWidgets.QLabel("", self.statusBar)
        self.statusLevel = QtWidgets.QLabel("", self.statusBar)
      
        self.statusBar.addWidget(self.statusPixmap, 0)
        self.statusBar.addWidget(self.__statusLabel, 0)
        self.statusBar.addWidget(self.progressbar, 1)
        self.statusBar.addWidget(self.statusSize, 0)
        self.statusBar.addWidget(self.statusZoom, 0)
        self.statusBar.addWidget(self.statusScale, 0)
        self.statusBar.addWidget(self.statusLevel, 0)


    def setMessage(self, message='', duration=0, pixmap='', size=0, zoom=0, scale=0, level=0, colour=''):
        self.__statusTimer.stop()
        self.lastMessage = str(self.__statusLabel.text())

        if len(message) > 0:
            self.__statusLabel.setText(message)
            self.__statusLabel.show()

        if size > 0:
            self.statusSize.setText(str(self.parent.size[0]) + ' x ' + str(self.parent.size[1]))
        
        if zoom > 0:
            self.statusZoom.setText(str(int(self.parent.factor)) + '%')
        
        if scale > 0:
            self.statusScale.setText("F="+str(self.parent.sarscale))
        
        if level > 0:
            self.statusLevel.setText("S="+str(self.parent.scale))
        
        if duration > 0:
            self.__statusTimer.start(1000 * duration, True)

        if pixmap:
            self.pixmapLabel.setPixmap(pixmap)
        
        if colour == 'G':
            self.statusPixmap.setPixmap(self.pixmap_green)
        elif colour == 'Y':
            self.statusPixmap.setPixmap(self.pixmap_yellow)
        elif colour == 'R':
            self.statusPixmap.setPixmap(self.pixmap_red)
        
        QtWidgets.QApplication.processEvents()
            
    def resetMessage(self):
        self.__statusTimer.stop()
        if self.lastMessage:
            self.__statusLabel.setText(self.lastMessage)
        else:
            self.__statusLabel.setText('')

