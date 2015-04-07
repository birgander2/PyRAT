import pyrat
import logging
import numpy as np
from PyQt4 import QtCore, QtGui
from pyrat.load.tools import rrat, RatFile, Xml2Py
from pyrat.filter.tools import rebin
from pyrat.tools import ProgressBar


from ipdb import set_trace
stop = set_trace


class RatCovar(pyrat.ImportWorker):

    #gui = {'menu': 'File', 'entry': 'Open RAT file', 'before': 'Open external'}
    #para = {
    #    'filename': {'value': '', 'type': 'openfile', 'text': ''}
    #}

    def __init__(self, *args, **kwargs):
        super(RatCovar, self).__init__(*args, **kwargs)
        self.name = "RAT COVAR IMPORT"        

        if 'rat_block' not in self.__dict__:
            self.rat_block = [0]
        if 'subx' not in self.__dict__:
            self.subx = 3
        if 'suby' not in self.__dict__:
            self.suby = 3

    def reader(self, *args, **kwargs):
        fN = len(self.filenames)

        meta = {}
        if 'ppfiles' in self.__dict__:
            pp = Xml2Py(self.ppfiles[0])
            meta['sensor'] = 'DLR F-SAR'
            meta['band'] = band
            meta['prf'] = pp.prf
            meta['c0'] = pp.c0
            meta['rd'] = pp.rd
            meta['rs'] = pp.rsf/subx
            meta['lam'] = pp.__dict__['lam']
            meta['band'] = pp.band
            meta['antdir'] = pp.antdir
            meta['v0'] = pp.v0
            meta['h0'] = pp.h0
            meta['terrain'] = pp.terrain
            meta['bw'] = pp.cbw
            if (len(self.rat_block) == 4):
                meta['rd'] += self.rat_block[0] / meta['rs']
            
            meta['CH_pol'] = [' ']*(fN**2)
            pind = 0
            for f1 in self.ppfiles:
                pp1 = Xml2Py(f1)
                for f2 in self.ppfiles:
                    pp2 = Xml2Py(f2)
                    meta['CH_pol'][pind] = pp1.polarisation+pp2.polarisation
                    ppind += 1

        file = RatFile(self.filenames[0])
        if len(self.rat_block) != 4:
            imDim = file.dim[0:2][::-1]
        else:
            imDim = self.rat_block[2:4][::-1]
        imDim = np.asarray(imDim,dtype='i4')

        blenOut = 100
        blenIn = blenOut*self.suby
        nBlocks = int(np.ceil(imDim[0]/float(blenIn)))

        sub = np.asarray([self.suby,self.subx],dtype='i4')
        oDim = imDim / sub

        meta['nrg'] = oDim[0]
        meta['naz'] = oDim[1]
            
        P = ProgressBar('  ' + self.name, nBlocks)
        P.update(0)
        for n in range(nBlocks):
            rInd = np.minimum(np.asarray([n*blenIn,(n+1)*blenIn]),imDim[0])
            wInd = rInd/sub
            if wInd[0] >= wInd[1]:
                continue

            binInd = [(wInd[1]-wInd[0])*sub[0],oDim[1]*sub[1]]
            bDim = (wInd[1]-wInd[0],oDim[1])
            for u in range(fN):
                blk = [0,rInd[0],imDim[1],rInd[1]-rInd[0]]
                if len(self.rat_block) == 4:
                    blk[0] = self.rat_block[0]
                    blk[1] += self.rat_block[1]
                    blk[2] = self.rat_block[2]
                arr1 = rrat(self.filenames[u],block=blk)
                arr1 = arr1[0:binInd[0],0:binInd[1]]

                if (u == 0) and (n == 0):
                    covar = np.empty((fN,fN)+tuple(oDim),dtype=arr1.dtype)

                covar[u,u,wInd[0]:wInd[1],:] = rebin(np.abs(arr1)**2,bDim)

                for v in range(u+1,fN):
                    arr2 = rrat(self.filenames[v],block=blk)
                    
                    covar[u,v,wInd[0]:wInd[1],:] = rebin(arr1*np.conj(arr2),bDim)
                    covar[v,u,wInd[0]:wInd[1],:] = np.conj(covar[u,v,wInd[0]:wInd[1],:])

            P.update(n+1)

        return covar, meta

    @classmethod
    def guirun_disabled(cls, viewer):
        filename = str(QtGui.QFileDialog().getOpenFileName())
        if filename:
            plugin = cls(filename=filename)
            plugin.run()
            del plugin
            viewer.statusBar.setMessage(message='Ready', colour='G')
            viewer.updateViewer()
