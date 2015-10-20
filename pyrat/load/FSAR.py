import pyrat
import glob, os
import logging
import copy
import numpy as np
from PyQt4 import QtGui, QtCore
from pyrat.load.tools import RatFile, Xml2Py
from pyrat.viewer.Dialogs import FlexFilesel
from pyrat.viewer.Widgets import HLine, CropBoxWidget, FileselWidget, ProductContentWidget


class FSAR(pyrat.ImportWorker):
    """
    Import of DLR F-SAR SLC product

    :param dir: The F-SAR product directory.
    :type dir: str
    :param match: A matching string to select subset of files
    :type match: string
    :param crop: A crop region / subset to import (az_start, az_end, rg_start, rg_end)
    :type crop: tuple
    :author: Andreas Reigber
    """
    gui = {'menu': 'File|Import airborne', 'entry': 'F-SAR'}
    para = [
        {'var': 'dir', 'value': ''},
        {'var': 'bands', 'value': '*'},
        {'var': 'polarisations', 'value': '*'},
        {'var': 'product', 'value': 'RGI-SLC'},
        {'var': 'crop', 'value': [0, 0, 0, 0]},
        {'var': 'mask', 'type': bool, 'value': False}]

    def __init__(self, *args, **kwargs):
        super(FSAR, self).__init__(*args, **kwargs)
        self.name = "FSAR SLC IMPORT"
        if 'crop' not in self.__dict__:
            self.crop = (0, 0, 0, 0)

    def reader(self, *args, **kwargs):

        if self.product == 'RGI-SLC':
            head = 'slc'
            src = ('RGI','RGI-SR')
        if self.product == 'RGI-AMP':
            head = 'amp'
            src = ('RGI','RGI-SR')
        if self.product == 'INF-SLC':
            head = 'slc_coreg'
            src = ('INF','INF-SR')
        if self.product == 'INF-CIR':
            head = 'slcpol'
            src = ('INF','INF-SR')

        if self.polarisations == '*':
            self.polarisations = '??'
        if self.product == 'INF-CIR':
            files = glob.glob(os.path.join(self.dir, src[0], src[1], head+'*'+self.bands.upper()+self.polarisations.lower()+'*_c'+str(self.track)+'_s'+str(self.subaperture)+'_*coreg.rat'))
        else:
            files = glob.glob(os.path.join(self.dir, src[0], src[1], head+'*'+self.bands.upper()+self.polarisations.lower()+'_*.rat'))
        bands = list(set([os.path.basename(slc).split('_')[2][0] for slc in files]))
        pols = list(set([os.path.basename(slc).split('_')[2][1:3] for slc in files]))

        # match = self.match if hasattr(self, 'match') else ''
        # files = glob.glob(self.filename+'/RGI/RGI-SR/slc_*'+match+'_*rat')
        # bands = set([f.split('_')[-2][0] for f in files])

        array = []
        meta = []
        for band in bands:
            if hasattr(self, 'band') and band not in self.band:
                logging.warning(band+'-band data not found in specified directory')
            else:
                bandfiles = [f for f in files if '_'+band in f]
                fil = RatFile(bandfiles[0])
                if self.mask is True:
                    maskfile = glob.glob(os.path.join(self.dir, src[0], src[1], 'mask*'+band.upper()+'*.rat'))
                    msk = RatFile(maskfile[0])
                naz = fil.dim[1]
                nrg = fil.dim[0]
                block = list(self.crop)
                if block[1] == 0:
                    block[1] = naz
                if block[3] == 0:
                    block[3] = nrg
                daz = block[1]-block[0]
                drg = block[3]-block[2]

                barr = np.empty((len(bandfiles), daz, drg), dtype='complex64')
                bmeta = {}
                bmeta['sensor'] = 'DLR F-SAR'
                bmeta['band'] = band
                bmeta['CH_pol'] = [' ']*len(bandfiles)
                for k, f in enumerate(bandfiles):
                    logging.info("Found "+f)
                    barr[k, ...] = RatFile(f).read(block=(block[2], block[0], drg, daz))
                    if self.mask is True:
                        mask = msk.read(block=(block[2], block[0], drg, daz))
                        # barr[k, ...] *= mask
                    if self.product == 'RGI-SLC':
                        ppfile = f.replace('RGI-SR', 'RGI-RDP').replace('slc_', 'pp_').replace('.rat', '.xml')
                    if self.product == 'RGI-AMP':
                        ppfile = f.replace('RGI-SR', 'RGI-RDP').replace('amp_', 'pp_').replace('.rat', '.xml')
                    if self.product == 'INF-SLC':
                        ppname = 'pp_'+'_'.join(os.path.basename(f).split('_')[3:]).replace('.rat','.xml')
                        ppfile = os.path.join(self.dir,'INF','INF-RDP',ppname)
                    if self.product == 'INF-CIR':
                        ppname = 'ppgeo_csar_'+'_'.join(os.path.basename(f).split('_')[1:4])+'.xml'
                        ppfile = os.path.join(self.dir,'GTC','GTC-RDP',ppname)
                    pp = Xml2Py(ppfile)
                    bmeta['CH_pol'][k] = pp.polarisation
                bmeta['prf'] = pp.prf
                bmeta['c0'] = pp.c0
                bmeta['rd'] = pp.rd
                bmeta['rsf'] = pp.rsf
                bmeta['nrg'] = drg
                bmeta['naz'] = daz
                bmeta['lam'] = pp.__dict__['lambda']
                bmeta['band'] = pp.band
                bmeta['antdir'] = pp.antdir
                bmeta['v0'] = pp.v0
                bmeta['bw'] = pp.cbw
                bmeta['ps_rg'] = pp.ps_rg
                bmeta['ps_az'] = pp.ps_az
                bmeta['rd'] += block[2] / bmeta['rsf']
                bmeta['h0'] = pp.h0
                bmeta['pre_az'] = pp.pre_az

                if self.mask is True:
                    array.append([barr, mask])
                else:
                    array.append(barr)
                meta.append(bmeta)


        if len(array) == 0:
            return None, None
        elif len(array) == 1:
            return array[0], meta[0]
        else:
            return array, meta

    @classmethod
    def guirun(cls, viewer):
        para_backup = copy.deepcopy(cls.para)                # keep a deep copy of the default parameters
        wid = FsarImportWidget()
        wid.update()
        res = wid.exec_()
        if res == 1:
            plugin = cls(dir=wid.dir, product=wid.product, bands=wid.bands, polarisations=wid.polar, crop=wid.crop)
            viewer.statusBar.setMessage(message=plugin.name+' running', colour = 'R')
            plugin.run()
            del plugin
            viewer.statusBar.setMessage(message='Ready', colour='G')
            viewer.updateViewer()


class FsarImportWidget(QtGui.QDialog):
    def __init__(self, parent=None, dir=None):
        super(FsarImportWidget, self).__init__(parent)
        self.setWindowTitle("FSAR import")
        mainlayout = QtGui.QVBoxLayout(self)

        self.dirwidget = FileselWidget(title='FSAR product dir', type='opendir')
        self.dirwidget.setvalue(dir)
        mainlayout.addWidget(self.dirwidget)
        mainlayout.addWidget(HLine())
        self.productwidget = ProductContentWidget(products=["RGI-SLC", "RGI-AMP", "INF-SLC"])
        mainlayout.addWidget(self.productwidget)
        mainlayout.addWidget(HLine())
        self.cropwidget = CropBoxWidget(title='Select crop (0=maximum)')
        mainlayout.addWidget(self.cropwidget)

        self.buttons = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
                                              QtCore.Qt.Horizontal, self)
        mainlayout.addWidget(self.buttons)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

        self.dirwidget.text.textChanged.connect(lambda: self.update(mode=0))   # update all

        self.bandupdate = lambda: self.update(mode=2)                          # update band
        self.productwidget.product.currentIndexChanged.connect(lambda: self.update(mode=1))# update product
        self.productwidget.band.currentIndexChanged.connect(self.bandupdate)

    def update(self, mode=0):
        self.dir = str(self.dirwidget.getvalue())
        self.product = self.productwidget.getvalue(0)
        self.bands = self.productwidget.getvalue(1)
        self.polar = self.productwidget.getvalue(2)

        if self.product == 'RGI-SLC':
            head = 'slc'
            src = ('RGI','RGI-SR')
            code_pos = 2
        if self.product == 'RGI-AMP':
            head = 'amp'
            src = ('RGI','RGI-SR')
            code_pos = 2
        if self.product == 'INF-SLC':
            head = 'slc_coreg'
            src = ('INF','INF-SR')
            code_pos = 4

        files = glob.glob(os.path.join(self.dir, src[0], src[1], head+'*'+self.bands.upper()+self.polar.lower()+'*.rat'))

        if mode == 0:
            allfiles = glob.glob(os.path.join(self.dir,src[0],src[1],head+'*.rat'))
            self.bands = '*'
            self.polar = '*'
        if mode == 1:
            allfiles = glob.glob(os.path.join(self.dir,src[0],src[1],head+'*.rat'))
        if mode == 2:
            allfiles = glob.glob(os.path.join(self.dir,src[0],src[1],head+'*'+self.bands.upper()+'*.rat'))

        # allfiles = glob.glob(os.path.join(self.dir,'RGI','RGI-SR',head+'*.rat'))
        allbands = list(set([os.path.basename(slc).split('_')[code_pos][0] for slc in allfiles]))
        allpols = list(set([os.path.basename(slc).split('_')[code_pos][1:3].upper() for slc in allfiles]))

        nrg = 0
        naz = 0
        for filename in files:
            lun = RatFile(filename)
            nrg = max(nrg, lun.dim[0])
            naz = max(naz, lun.dim[1])
        self.cropwidget.setrange([[0, naz], [0, naz], [0, nrg], [0, nrg]])
        self.cropwidget.setvalues([0, naz, 0, nrg])

        if mode == 0 or mode == 1:
            self.productwidget.band.currentIndexChanged.disconnect(self.bandupdate)
            self.productwidget.updatepolar(allpols)
            self.productwidget.updatebands(allbands)
            self.productwidget.setvalue(1, self.bands)
            self.productwidget.setvalue(2, self.polar)
            self.productwidget.band.currentIndexChanged.connect(self.bandupdate)

        elif mode == 2:
            self.productwidget.updatepolar(allpols)
            self.productwidget.setvalue(2, self.polar)

    def accept(self):
        self.product = self.productwidget.getvalue(0)
        self.bands = self.productwidget.getvalue(1)
        self.polar = self.productwidget.getvalue(2)
        self.crop = self.cropwidget.getvalues()
        super(FsarImportWidget, self).accept()


def fsar(*args, **kwargs):
    return FSAR(*args, **kwargs).run(**kwargs)

