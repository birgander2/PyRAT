from yapsy.PluginManager import PluginManager
from PyQt4 import QtGui, QtCore
from pylab import *
import PyRat, scipy, numpy
import IDL
from StatusBar import *
import pdb

class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        self.box      = [0,100,0,100]
        self.size     = [100,100]
        self.factor   = 1.0
        self.sarscale = 2.5
        self.type     = 'A'

        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle("PyRAT - Radar Tools")
        self.makeActions()
        self.makeToolbar()
        self.makeStatusBar()
        self.makeMenu()
        self.makeView()
        self.initPlugins()
        self.resize(800, 800)
        self.show()

#-------------------------------- TOOL BAR
    def makeToolbar(self):
        self.openTB  = QtGui.QAction(QtGui.QIcon('icons/document-open.png'), 'Open', self)
        self.closeTB = QtGui.QAction(QtGui.QIcon('icons/document-close.png'), 'Close', self)
        #self.zoominTB = QtGui.QAction(QtGui.QIcon('icons/zoom-in.png'), 'Zoom in', self)
        #self.zoomoutTB = QtGui.QAction(QtGui.QIcon('icons/zoom-out.png'), 'Zoom out', self)
        #self.zoomresetTB = QtGui.QAction(QtGui.QIcon('icons/zoom-fit-best.png'), 'Fit zoom', self)
        self.seperatorTB = QtGui.QAction(self)

        self.toolbar1 = self.addToolBar("File")
        self.toolbar1.addAction(self.openTB)       
        self.toolbar1.addAction(self.closeTB)    

        self.toolbar2 = self.addToolBar("Display")
        self.toolbar2.addAction(self.zoomOutAct)       
        self.viewCombo = QtGui.QComboBox(self)
        self.viewCombo.insertItems(1, ["100%", "Fit to window","Fit to width","Fit to height", "100%"])
        self.viewCombo.setEditable(False)
        self.viewCombo.activated.connect(self.comboZoom)
        self.toolbar2.addWidget(self.viewCombo)
        self.toolbar2.addAction(self.zoomInAct)       
       
        #self.toolbar3 = self.addToolBar("Layer")
        #self.toolbar3.addAction(self.zoominTB)       
        #self.toolbar3.addAction(self.zoomresetTB)       
       
#-------------------------------- STATUS BAR
    def makeStatusBar(self):   
        self.statusBar = StatusBar(self)
        self.statusBar.setMessage(message='no data loaded')
        
#-------------------------------- DEFAULT ACTIONS (only those not implemented trough plugins)
    def makeActions(self):
        self.exitAct        = QtGui.QAction('Exit', self, shortcut='Q', triggered = self.close)       
        self.zoomInAct      = QtGui.QAction(QtGui.QIcon('icons/zoom-in.png'), "Zoom &In (25%)", self, shortcut="up", triggered=lambda: self.zoom(3.0/2.0))
        self.zoomOutAct     = QtGui.QAction(QtGui.QIcon('icons/zoom-out.png'),"Zoom &Out (25%)", self, shortcut="down", triggered=lambda: self.zoom(2.0/3.0))
        self.fitToWindowAct = QtGui.QAction(QtGui.QIcon('icons/zoom-fit-best.png'), "Reset view", self, shortcut="f", triggered=self.resetView)
        self.viewAmpAct     = QtGui.QAction("View as amplitude", self, checkable=True, shortcut="1", triggered=self.viewAsAmplitude)
        self.viewPhaAct     = QtGui.QAction("View as phase", self, checkable=True, shortcut="2", triggered=self.viewAsPhase)
        self.viewCohAct     = QtGui.QAction("View as coherence", self, checkable=True, shortcut="3", triggered=self.viewAsCoherence)
        self.viewBrighter   = QtGui.QAction("View brighter", self, shortcut="right", triggered=self.brighterView)
        self.viewDarker     = QtGui.QAction("View darker", self, shortcut="left", triggered=self.darkerView)
 
#---------------------------------- PLUGINS
    def initPlugins(self):
        manager = PluginManager()
        manager.setPluginPlaces(["plugins","PyRat/Import","PyRat/Help"])
        manager.collectPlugins()
        print "-> Registering plugins into GUI"
        for self.plugin in manager.getAllPlugins():    # still need to evaluate plugin.order... 
            print self.plugin.name, self.plugin.menu, self.plugin.order
            self.plugin.plugin_object.registerGUI(self)

#---------------------------------- MENUBAR
    def makeMenu(self):
        self.menubar = self.menuBar()
        self.menue = {
        "File"  : self.menubar.addMenu('File'),
        "View"  : self.menubar.addMenu('View'), 
        "Tools" : self.menubar.addMenu('Tools'),
        "SAR"   : self.menubar.addMenu('SAR'),
        "PolSAR": self.menubar.addMenu('PolSAR'),
        "Help"  : self.menubar.addMenu('Help')}
        
        self.menue["File"].addSeparator()
        self.menue["File"].addAction(self.exitAct)
        
        self.menue["View"].addAction(self.fitToWindowAct)
        self.menue["View"].addAction(self.zoomInAct)
        self.menue["View"].addAction(self.zoomOutAct)
        self.menue["View"].addSeparator()
        self.menue["View"].addAction(self.viewBrighter)
        self.menue["View"].addAction(self.viewDarker)
        self.menue["View"].addSeparator()
        self.viewSel = QtGui.QActionGroup(self.menue["View"],exclusive=True)
        foo = self.viewSel.addAction(self.viewAmpAct)
        self.menue["View"].addAction(foo)
        foo = self.viewSel.addAction(self.viewPhaAct)
        self.menue["View"].addAction(foo)
        foo = self.viewSel.addAction(self.viewCohAct)
        self.menue["View"].addAction(foo)
        
        #self.menue["View"].addAction(self.viewAmpAct)
        #self.menue["View"].addAction(self.viewPhaAct)
        #self.menue["View"].addAction(self.viewCohAct)
     
#---------------------------------- VIEW AREA
    def makeView(self):
        self.imageLabel = QtGui.QLabel()
        self.imageLabel.setBackgroundRole(QtGui.QPalette.Base)
        self.imageLabel.setStyleSheet("QLabel { background-color: #333 }")
        self.imageLabel.setSizePolicy(QtGui.QSizePolicy.Ignored,QtGui.QSizePolicy.Ignored)
        self.imageLabel.setScaledContents(False)
        self.setCentralWidget(self.imageLabel)

#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
#  CODE FROM RATVIEWER
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
#------------------------------------DISPLAY DATA
    def genPyramid(self, data):
        if data.ndim > 2:        # Temporay workaround
            data = data[0,...]   # for multichannel data
        
        if self.type == 'P' and numpy.iscomplexobj(data):
            data = [numpy.angle(data)]
        elif self.type == 'A' and numpy.iscomplexobj(data):
            data = [numpy.abs(data)]
        else:
            data = [data]
        data[0][isnan(data[0])] = 0.0
        scale = 0
        lange = min(data[scale].shape)
        while lange > 1:
            #self.statusBar.setMessage('Generating pyramid layer : '+str(scale+1))
            cut = numpy.array(data[scale].shape)//2*2
            if self.type == 'P': 
                data.append(IDL.rebin(data[scale][:cut[0],:cut[1]],numpy.array(data[scale].shape)//2,phase=True))
            else: 
                data.append(IDL.rebin(data[scale][:cut[0],:cut[1]],numpy.array(data[scale].shape)//2))

            scale += 1
            lange = min(data[scale].shape)
        self.data = data
        self.size = data[0].shape[::-1]   # [xmin,xmax,ymin,ymax]
    
    def data2img(self, cut_box, scale=0):
        if self.type == 'P':   
            img = uint8(clip(self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1]]/pi*128+127,0.0,255.0))
            self.viewPhaAct.setChecked(True)
        elif self.type == 'C': 
            img = uint8(clip(self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1]]*255.0,0.0,255.0))
            self.viewCohAct.setChecked(True)
        elif self.type == 'A': 
            img = sarscl(self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1]],factor=self.sarscale)
            self.viewAmpAct.setChecked(True)
        else:
            img = uint8(self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1]])

        img = img[0:img.shape[0]//4*4,0:img.shape[1]//4*4]         # QT Limitation, better move to cut_box calculation
        return QtGui.QImage(img.tostring(), img.shape[1], img.shape[0], QtGui.QImage.Format_Indexed8)   # convert to Qimage
        
    def processPicture(self, **kwargs):
        
        if "fitwin" in kwargs.keys():                        
            self.scale = len(self.data)-1
            while self.scale >= 0:
                if  self.data[self.scale].shape[1] > self.imageLabel.width() or self.data[self.scale].shape[0] > self.imageLabel.height(): break
                self.scale -= 1
            self.box = [0,self.size[0],0,self.size[1]]
            cut_box = [foo // 2**self.scale for foo in self.box]
        else:
            self.scale = len(self.data)-1    
            while self.scale >= 0:
                cut_box = [foo // 2**self.scale for foo in self.box]
                #print cut_box
                if  cut_box[1]-cut_box[0] > self.imageLabel.width() or cut_box[3]-cut_box[2] > self.imageLabel.height(): break
                self.scale -= 1
        self.scale = clip(self.scale,0,len(self.data)-1)       
         
        img = self.data2img(cut_box, scale=self.scale)
        
        xWin = self.imageLabel.width()
        yWin = self.imageLabel.height()
        winRatio = 1.0*xWin/yWin

        self.width  = img.width()
        self.height = img.height()
        imgRatio = 1.0*self.width/self.height

        if imgRatio >= winRatio: #match widths
                self.width  = xWin
                self.height = xWin/imgRatio
        else:   #match heights
                self.height = yWin
                self.width  = yWin*imgRatio

        self.factor = int(100.0*self.width/(self.box[1]-self.box[0]))
        if self.factor <= 100:
            img = img.scaled(int(self.width),int(self.height))  # Bilinear?
        else:
            img = img.scaled(int(self.width),int(self.height))  # Nearest Neighbour
        
        
        self.statusBar.setMessage(size=1, zoom=1, level=1, scale=1)
        self.viewCombo.setItemText(0,str(int(self.factor)) + '%')
        
        colortable = [QtGui.qRgb(i,i,i) for i in xrange(256)]
        img.setColorTable(colortable)
        self.imageLabel.setPixmap(QtGui.QPixmap.fromImage(img))
    
    def zoom(self,factor,mx=0,my=0):
        px = self.box[0] + int((1.0*self.box[1]-self.box[0])/self.imageLabel.width()*(mx+self.imageLabel.width()//2))
        py = self.box[2] + int((1.0*self.box[3]-self.box[2])/self.imageLabel.height()*(my+self.imageLabel.height()//2))
        
        midx = self.box[0]+(self.box[1]-self.box[0])//2 
        midy = self.box[2]+(self.box[3]-self.box[2])//2 
        sizx = self.box[1]-self.box[0]
        sizy = self.box[3]-self.box[2]
        newx = clip(int(sizx/factor),0,self.size[0])
        newy = clip(int(sizy/factor),0,self.size[1])
        imgRatio = 1.0*newx/newy
        xWin = self.imageLabel.width()
        yWin = self.imageLabel.height()
        winRatio = 1.0*xWin/yWin
        if imgRatio >= winRatio: #match widths:
            newy = clip(int(newx / winRatio),0,self.size[1])
        else:  
            newx = clip(int(newy * winRatio),0,self.size[0])
        newx = clip(newx,4,self.size[0])
        newy = clip(newy,4,self.size[1])
        if midx-newx//2 < 0:            midx = newx//2       
        if midx+newx//2 > self.size[0]: midx = self.size[0]-newx//2       
        if midy-newy//2 < 0:            midy = newy//2       
        if midy+newy//2 > self.size[1]: midy = self.size[1]-newy//2       
        self.box = [midx-newx//2,midx+newx//2,midy-newy//2,midy+newy//2]

        if mx != 0 and my != 0: 
            midx = px - int((1.0*self.box[1]-self.box[0])/self.imageLabel.width()*mx)
            midy = py - int((1.0*self.box[3]-self.box[2])/self.imageLabel.height()*my)
            sizx = self.box[1]-self.box[0]
            sizy = self.box[3]-self.box[2]
            if midx-sizx//2 < 0:            midx = sizx//2       
            if midx+sizx//2 > self.size[0]: midx = self.size[0]-sizx//2       
            if midy-sizy//2 < 0:            midy = sizy//2       
            if midy+sizy//2 > self.size[1]: midy = self.size[1]-sizy//2       
            self.box = [midx-sizx//2,midx+sizx//2,midy-sizy//2,midy+sizy//2]
        if hasattr(self,'data'): self.processPicture()

    def resizeEvent(self, event):
        midx = self.box[0]+(self.box[1]-self.box[0])//2
        midy = self.box[2]+(self.box[3]-self.box[2])//2
        sizx = self.box[1]-self.box[0]
        sizy = self.box[3]-self.box[2]
        newx = clip(int(sizx),0,self.size[0])
        newy = clip(int(sizy),0,self.size[1])
        imgRatio = 1.0*newx/newy
        xWin = self.imageLabel.width()
        yWin = self.imageLabel.height()
        winRatio = 1.0*xWin/yWin
        if imgRatio >= winRatio: #match widths:
            newy = clip(int(newx / winRatio),0,self.size[1])
        else:  
            newx = clip(int(newy * winRatio),0,self.size[0])  
        self.box = [midx-newx//2,midx+newx//2,midy-newy//2,midy+newy//2]
        if hasattr(self,'data'): self.processPicture()
                
    def wheelEvent(self,event):
        if event.delta() < 0: self.zoom(2.0/3.0,mx=event.x()-self.imageLabel.width()/2,my=event.y()-self.imageLabel.height()/2)       
        if event.delta() > 0: self.zoom(3.0/2.0,mx=event.x()-self.imageLabel.width()/2,my=event.y()-self.imageLabel.height()/2)       
        if hasattr(self,'data'): self.processPicture()
       
    def mousePressEvent(self,event):
        self.dragX = event.x()      
        self.dragY = event.y() 
        
    def mouseMoveEvent(self, event):   
        pass
    
    def viewAsAmplitude(self):
        self.type = 'A'
        if hasattr(self,'data'): self.processPicture()

    def viewAsCoherence(self):
        self.type = 'C'
        if hasattr(self,'data'): self.processPicture()

    def viewAsPhase(self):
        self.type = 'P'
        if hasattr(self,'data'): self.processPicture()

    def darkerView(self):
        self.sarscale += 0.5
        if hasattr(self,'data'): self.processPicture()

    def brighterView(self):
        self.sarscale -= 0.5
        if self.sarscale < 0.5: self.sarscale = 0.5
        if hasattr(self,'data'): self.processPicture()

    def resetView(self):
        self.sarscale = 2.5
        if hasattr(self,'data'): self.processPicture(fitwin=1)

    def comboZoom(self, index):
        if index == 1:
            if hasattr(self,'data'): self.processPicture(fitwin=1)
            self.viewCombo.setCurrentIndex(0)
        elif index == 2:
            print "Not implemented"
            self.viewCombo.setCurrentIndex(0)
        elif index == 3:
            print "Not implemented"
            self.viewCombo.setCurrentIndex(0)
        elif index == 4:
            print "Not implemented"
            self.viewCombo.setCurrentIndex(0)

    def mouseReleaseEvent(self,event):
        if event.button() == 1:
            dx = int((1.0*self.box[1]-self.box[0])/self.imageLabel.width()*(self.dragX - event.x()))
            dy = int((1.0*self.box[3]-self.box[2])/self.imageLabel.height()*(self.dragY - event.y()))
            midx = self.box[0]+(self.box[1]-self.box[0])//2 + dx
            midy = self.box[2]+(self.box[3]-self.box[2])//2 + dy
            sizx = self.box[1]-self.box[0]
            sizy = self.box[3]-self.box[2]
            if midx-sizx//2 < 0:            midx = sizx//2       
            if midx+sizx//2 > self.size[0]: midx = self.size[0]-sizx//2       
            if midy-sizy//2 < 0:            midy = sizy//2       
            if midy+sizy//2 > self.size[1]: midy = self.size[1]-sizy//2       
            self.box = [midx-sizx//2,midx+sizx//2,midy-sizy//2,midy+sizy//2]
            if hasattr(self,'data'): self.processPicture()
    
def sarscl(img, factor=2.5):
    return uint8(clip(255.0*img/factor/mean(img[img > 0],dtype='f8'),0,255))
   
 