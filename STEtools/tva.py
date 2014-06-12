try:
    from PyQt4 import QtGui, QtCore
    USE_PYSIDE = False
except ImportError:
    from PySide import QtCore, QtGui
    USE_PYSIDE = True
import STEtools as STE
import IDL
import numpy
import __builtin__

class SARImageView(QtGui.QGraphicsView):

    def __init__(self, parent=None, name="SARImageView", *args, **kargs):
        self.box      = [0,100,0,100]
        self.size     = [100,100]
        self.factor   = 1.0
        self.sarscale = 2.5
        self.type     = 'A'
        self.colour   = False
        if 'type' in kargs: self.type = kargs['type'].upper()
        if self.type not in ['A','C','P','B']: self.type     = 'A'

        if hasattr(__builtin__,'myPalette'):
            self.colortable = [QtGui.qRgb(__builtin__.myPalette[i,0],__builtin__.myPalette[i,1],__builtin__.myPalette[i,2]) for i in xrange(256)]
        else:
            self.colortable = [QtGui.qRgb(i,i,i) for i in xrange(256)]
        
        QtGui.QGraphicsView.__init__(self, parent)
        self.setGeometry(QtCore.QRect(0, 0, 520, 520))
        self.setObjectName("SARImageView")
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scene = QtGui.QGraphicsScene(self)
        self.setBackgroundBrush(QtGui.QColor(50,50,50))
        self.scene.setSceneRect(QtCore.QRectF(0, 0, 512, 512))
        self.setScene(self.scene)
                   
    def widget_close(self):
        """Closes the widget nicely, making sure to clear the graphics scene and release memory."""
        self.scene.clear()
        self.close()
        del self.data
        self.setParent(None)
        
    def genPyramid(self, data, **kargs):
        if data.ndim == 3:
            self.colour = True
            dshape = data.shape
            dmin = dshape.index(min(dshape))
            if dmin == 0:
                data = numpy.rollaxis(numpy.rollaxis(data,2),2)    
            elif dmin == 1:
                data = numpy.rollaxis(data,2,1)
            else:
                pass
            dshape = data.shape
            if dshape[2] < 3:
                data = numpy.append(data,data[...,0:1],axis=2)
            if dshape[2] > 3:
                data = data[...,0:3]

        if self.type == 'P' and numpy.iscomplexobj(data):
            data = [numpy.angle(data)]
        elif numpy.iscomplexobj(data):
            data = [numpy.abs(data)]
        else:
            data = [data]
        data[0][numpy.isnan(data[0])] = 0.0
        scale = 0
        lange = min(data[scale].shape)
        while lange > 1:
            cut = numpy.array(data[scale].shape[0:2])//2*2
            newdim = numpy.array(data[scale].shape[0:2])//2
            if self.colour == True:
                newdim = numpy.append(newdim,3)
            if self.type == 'P': 
                data.append(IDL.rebin(data[scale][:cut[0],:cut[1],...],newdim,phase=True))
            else: 
                data.append(IDL.rebin(data[scale][:cut[0],:cut[1],...],newdim))
            scale += 1
            lange = min(data[scale].shape)
        self.data = data
        self.size = data[0].shape[0:2][::-1]   # [xmin,xmax,ymin,ymax]
            
    def data2img(self, cut_box, scale=0):        
        if self.data[scale].dtype == 'uint8':
            img = self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1],...]
        elif self.type == 'P':
            img = STE.phascale(self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1],...])
        elif self.type == 'C': 
            img = STE.cohscale(self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1],...])
        elif self.type == 'A':
            img = self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1],...].copy()
            if self.colour == True:
                for k in range(img.shape[2]): img[...,k] /= numpy.mean(img[...,k])
            img = STE.sarscale(img,factor=self.sarscale)
        elif self.type == 'B': 
            img = numpy.uint8(numpy.clip(self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1],...],0,255))
        else:
            img = IDL.bytscl(self.data[scale][cut_box[2]:cut_box[3],cut_box[0]:cut_box[1],...])
        img = img[0:img.shape[0]//4*4,0:img.shape[1]//4*4,...]         # QT Limitation, better move to cut_box calculation
        
        if self.colour == True:
            return QtGui.QImage(img.tostring(), img.shape[1], img.shape[0], QtGui.QImage.Format_RGB888)   # convert to Qimage
        else:
            return QtGui.QImage(img.tostring(), img.shape[1], img.shape[0], QtGui.QImage.Format_Indexed8)   # convert to Qimage
        
    def processPicture(self, *args, **kwargs):
        
        if "fitwin" in kwargs.keys():                        
            self.scale = len(self.data)-1
            while self.scale >= 0:
                if  self.data[self.scale].shape[1] > self.rect().width() or self.data[self.scale].shape[0] > self.rect().height(): break
                self.scale -= 1
            self.box = [0,self.size[0],0,self.size[1]]
            cut_box = [foo // 2**self.scale for foo in self.box]
        else:
            self.scale = len(self.data)-1    
            while self.scale >= 0:
                cut_box = [foo // 2**self.scale for foo in self.box]
                if  cut_box[1]-cut_box[0] > self.rect().width() or cut_box[3]-cut_box[2] > self.rect().height(): break
                self.scale -= 1
        self.scale = numpy.clip(self.scale,0,len(self.data)-1)  
        img = self.data2img(cut_box, scale=self.scale)
        
        xWin = self.rect().width()
        yWin = self.rect().height()
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
        
        img.setColorTable(self.colortable)
        
        if hasattr(self,'imageItem'): self.scene.removeItem(self.imageItem)
        self.scene.setSceneRect(QtCore.QRectF(img.rect()))
        self.imageItem = self.scene.addPixmap(QtGui.QPixmap.fromImage(img))
        
        #self.setDragMode(QtGui.QGraphicsView.ScrollHandDrag)
        #pdb.set_trace()
        
    def zoom(self,factor,mx=0,my=0):
        px = self.box[0] + int((1.0*self.box[1]-self.box[0])/self.rect().width()*(mx+self.rect().width()//2))
        py = self.box[2] + int((1.0*self.box[3]-self.box[2])/self.rect().height()*(my+self.rect().height()//2))
        
        midx = self.box[0]+(self.box[1]-self.box[0])//2 
        midy = self.box[2]+(self.box[3]-self.box[2])//2 
        sizx = self.box[1]-self.box[0]
        sizy = self.box[3]-self.box[2]
        newx = numpy.clip(int(sizx/factor),0,self.size[0])
        newy = numpy.clip(int(sizy/factor),0,self.size[1])
        imgRatio = 1.0*newx/newy
        xWin = self.rect().width()
        yWin = self.rect().height()
        winRatio = 1.0*xWin/yWin
        if imgRatio >= winRatio:  # match widths:
            newy = numpy.clip(int(newx / winRatio),0,self.size[1])
        else:  
            newx = numpy.clip(int(newy * winRatio),0,self.size[0])
        newx = numpy.clip(newx,4,self.size[0])
        newy = numpy.clip(newy,4,self.size[1])
        if midx-newx//2 < 0:            midx = newx//2       
        if midx+newx//2 > self.size[0]: midx = self.size[0]-newx//2       
        if midy-newy//2 < 0:            midy = newy//2       
        if midy+newy//2 > self.size[1]: midy = self.size[1]-newy//2       
        self.box = [midx-newx//2,midx+newx//2,midy-newy//2,midy+newy//2]

        if mx != 0 and my != 0: 
            midx = px - int((1.0*self.box[1]-self.box[0])/self.rect().width()*mx)
            midy = py - int((1.0*self.box[3]-self.box[2])/self.rect().height()*my)
            sizx = self.box[1]-self.box[0]
            sizy = self.box[3]-self.box[2]
            if midx-sizx//2 < 0:            midx = sizx//2       
            if midx+sizx//2 > self.size[0]: midx = self.size[0]-sizx//2       
            if midy-sizy//2 < 0:            midy = sizy//2       
            if midy+sizy//2 > self.size[1]: midy = self.size[1]-sizy//2       
            self.box = [midx-sizx//2,midx+sizx//2,midy-sizy//2,midy+sizy//2]
        if hasattr(self,'data'): self.processPicture()

    def savePicture(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save window content as JPG', '.')
        if filename:
            self.imageItem.pixmap().save(filename)
        
    def saveFullPicture(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save entire data set as JPG', '.')
        if filename:
            img = self.data2img([0,self.size[0],0,self.size[1]])
            img.setColorTable(self.colortable)
            img.save(filename)
            
    def resizeEvent(self, event):
        midx = self.box[0]+(self.box[1]-self.box[0])//2
        midy = self.box[2]+(self.box[3]-self.box[2])//2
        sizx = self.box[1]-self.box[0]
        sizy = self.box[3]-self.box[2]
        newx = numpy.clip(int(sizx),0,self.size[0])
        newy = numpy.clip(int(sizy),0,self.size[1])
        imgRatio = 1.0*newx/newy
        xWin = self.rect().width()
        yWin = self.rect().height()
        winRatio = 1.0*xWin/yWin
        if imgRatio >= winRatio: #match widths:
            newy = numpy.clip(int(newx / winRatio),0,self.size[1])
        else:  
            newx = numpy.clip(int(newy * winRatio),0,self.size[0])  
        self.box = [midx-newx//2,midx+newx//2,midy-newy//2,midy+newy//2]
        if hasattr(self,'data'): self.processPicture()

    def wheelEvent(self,event):
        if event.delta() < 0: self.zoom(2.0/3.0,mx=event.x()-self.rect().width()/2,my=event.y()-self.rect().height()/2)       
        if event.delta() > 0: self.zoom(3.0/2.0,mx=event.x()-self.rect().width()/2,my=event.y()-self.rect().height()/2)       
        if hasattr(self,'data'): self.processPicture()
       
    def mousePressEvent(self,event):
        self.dragX = event.x()      
        self.dragY = event.y()
        if event.button() == 2:
            posx = self.box[0]+int((1.0*self.box[1]-self.box[0])/self.rect().width()*self.dragX)
            posy = self.box[2]+int((1.0*self.box[3]-self.box[2])/self.rect().height()*self.dragY)
            if self.type == 'P': 
                val  = str(180.0/numpy.pi*self.data[0][posy,posx])+'deg'
            else: 
                val  = str(self.data[0][posy,posx])
            self.textItem =  self.scene.addText('')
            self.textItem.setHtml(str("<div style='background-color: #ffff00;'>") + 'x='+str(posx)+' y='+str(posy)+' val='+val+"</div>");
        
    def mouseReleaseEvent(self,event):
        if event.button() == 1:
            dx = int((1.0*self.box[1]-self.box[0])/self.rect().width()*(self.dragX - event.x()))
            dy = int((1.0*self.box[3]-self.box[2])/self.rect().height()*(self.dragY - event.y()))
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
        if event.button() == 2:
            if hasattr(self,'textItem'): self.scene.removeItem(self.textItem)
 
    def keyPressEvent(self,event):
        code = event.key()
        if code == QtCore.Qt.Key_Up: #up
            self.zoom(3.0/2.0)
            if hasattr(self,'data'): self.processPicture()
        elif code == QtCore.Qt.Key_Down:
            self.zoom(2.0/3.0)
            if hasattr(self,'data'): self.processPicture()
        elif code == QtCore.Qt.Key_Left:
            self.sarscale += 0.5
            if hasattr(self,'data'): self.processPicture()
        elif code == QtCore.Qt.Key_Right:
            self.sarscale -= 0.5
            if self.sarscale < 0.5: self.sarscale = 0.5
            if hasattr(self,'data'): self.processPicture()           
        elif code == QtCore.Qt.Key_F:  #F = FIT
            self.sarscale = 2.5
            if hasattr(self,'data'): self.processPicture(fitwin=1)
        elif code == QtCore.Qt.Key_1:  #1 = Amplitude
            self.type = 'A'
            if hasattr(self,'data'): self.processPicture()           
        elif code == QtCore.Qt.Key_2:  #2 = Phase
            self.type = 'P'
            if hasattr(self,'data'): self.processPicture()           
        elif code == QtCore.Qt.Key_3:  #3 = Coherence
            self.type = 'C'
            if hasattr(self,'data'): self.processPicture()           
        elif event.text() == 's':  #s = SAVE
            self.savePicture()
        elif event.text() == 'S':  #S = SAVE FULL RES
            self.saveFullPicture()
        
class SARImageWindow(SARImageView):
    def __init__(self, *args, **kargs):
        global QAPP
        inst = QtGui.QApplication.instance()
        if inst is None:
            QAPP = QtGui.QApplication([])
        else:
            QAPP = inst
        self.win = QtGui.QMainWindow()
        self.xwinsize = 512
        self.ywinsize = 512
        if 'xsize' in kargs: self.xwinsize = kargs['xsize']
        if 'ysize' in kargs: self.ywinsize = kargs['ysize']
        if 'size' in kargs:  self.ywinsize, self.xwinsize = kargs['size']
        self.win.resize(self.xwinsize,self.ywinsize)
        SARImageView.__init__(self, self.win, *args, **kargs)
        self.update(*args, **kargs)
        
    def update(self, *args, **kargs):
        if 'title' in kargs:
            self.win.setWindowTitle(kargs['title'])
            #del kargs['title']
        if len(args) > 0 or len(kargs) > 0:
            self.genPyramid(*args, **kargs)
            kargs['fitwin'] = True
            self.processPicture(*args, **kargs)
        self.win.setCentralWidget(self)
        self.win.show()
        QtCore.QCoreApplication.processEvents()   # needed to do screen updates!

def tva(*args, **kargs):
    """
    Opens a SARImageWindow with a SAR image visualised.
    
    Keyboard commands:  Zoom in/out, brighter/darker with the cursor keys
                        Switch between amplitude, phase and coherence display: 1,2,3
                        Save window / full image: s,S
    Mouse commands:     Drag with left button to move image
                        Scroll wheel to zoom in/out
                        Right button: Show current coordinate
                        
    Command line options:  xsize=900, ysize=200, size=(900,200)
                        Used for setting the initial window size, size 
                        overrides the other two.
    """
    if 'win' in kargs:
        kargs['win'].update(*args, **kargs)
        #del kargs['win']
    else:
        w = SARImageWindow(*args, **kargs)
        #w.show()
        return w

def tvp(*args, **kargs):
    """
    Opens a SARImageWindow with a phase image visualised
    """
    w = SARImageWindow(*args, type='P', **kargs)
    w.show()
    return w

def tvc(*args, **kargs):
    """
    Opens a SARImageWindow with a coherence image visualised
    """
    w = SARImageWindow(*args, type='C', **kargs)
    w.show()
    return w

