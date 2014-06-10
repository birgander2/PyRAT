import PyRat
import logging, copy, time
import numpy as np
import STEtools as STE
import multiprocessing as mp

def exec_parallel(args):
    return args[0].filter(args[1], meta=args[2], block=args[3]), args[2]
    
class FilterWorker(object):
    def __init__(self, *args, **kwargs):
        for (k, v) in kwargs.items():  # copy keywords to self
            setattr(self, k, v)
        self.name = 'UNKNOWN'
        self.blockprocess  = False
        self.allowed_dtype = False 
        self.allowed_ndim  = False
        self.require_para  = False
        self.blocksize     = 128
        self.blockoverlap  = 0
        self.nthreads      = PyRat.MP_Pool._processes
        self.noshow = ['noshow','nthreads','blockoverlap','blocksize', 'name','blockprocess','allowed_dtype','allowed_ndim','require_para']
    
    def filter(self, data, *args, **kwargs):
        logging.error(self.name + ': No filter method defined')
        return False
    
    def run_single(self, *args, **kwargs):
        print "SINGLE PROCESS - WARNING: NOT UP-TO-DATE"
        logging.info(self.name + '  '+str(dict((k, v) for k,v in self.__dict__.items() if k not in self.noshow)))
        if self.checkinput():
            meta  = PyRat.Data.getAnnotation()
            track = PyRat.Data.getTrack()
            result   = self.filter(PyRat.Data.getData(), track=track, meta=meta, *args, **kwargs)
            newlayer = PyRat.Data.addLayer(result, *args, **kwargs)
            PyRat.Data.setAnnotation(meta, layers=newlayer)
            PyRat.Data.setTrack(track, layer=newlayer)
            return newlayer
        else:
            return False
        
    def run(self, *args, **kwargs):
        logging.info(self.name + '  '+str(dict((k, v) for k,v in self.__dict__.items() if k not in self.noshow)))
        if self.checkinput():
            self.pre()
            if 'replace' in kwargs and self.blockprocess==True:
                pass
                #logging.warning("Blockprocessing and replace keyword not compatible - ignoring it!")
                #del kwargs['replace']
            self.il   = PyRat.Data.active
            self.ol   = None
            meta      = PyRat.Data.getAnnotation(layer=self.il)
            self.ydim = PyRat.Data.dshape[0]
            
            if isinstance(self.ydim, tuple):
                self.ydim = self.ydim[0]
            sidx = self.calc_blocks(self.ydim)

            if len(sidx) > 1 and self.nthreads > 1: 
                idx = [sidx[i:i+self.nthreads] for i in range(0, len(sidx), self.nthreads)]
            else:
                idx  = sidx

            if len(sidx) > 1:
                P = STE.ProgressBar('  '+self.name, len(sidx)-1)
                P.update(0)
                for k in range(len(idx)):
                    metain = copy.deepcopy(meta)
                    if self.nthreads > 1:
                        bidx = idx[k]
                        inputs = []
                        for l in range(len(bidx)):
                            block= (bidx[l],bidx[l]+self.blocksize,0,0)
                            data = self.read_block(l, bidx)
                            inputs.append((self, data, metain, block))
                        result = PyRat.MP_Pool.imap(exec_parallel, inputs)
                        for l, res in enumerate(result):
                            metain = res[1]
                            pos = sidx.index(bidx[l])
                            if res[0] != None:
                                self.save_block(res[0], pos, sidx, *args, **kwargs)
                            P.update(pos)
                    else:  # only 1 thread, but still blockprocessing
                        block= PyRat.Data.calcBlock((sidx[k],sidx[k]+self.blocksize,0,0))
                        data = self.read_block(k, sidx)
                        result = self.filter(data, meta=metain, block=block, *args, **kwargs)
                        if result != None:
                            self.save_block(result, k, sidx) 
                        if len(sidx) > 1: P.update(k)
            else: # no blockprocessing
                metain = copy.deepcopy(meta)
                block  = (0,self.blocksize,0,0)
                result  = self.filter(PyRat.Data.getData(), meta=metain, block=block, *args, **kwargs)
                if result != None:
                    self.ol = PyRat.Data.addLayer(result, *args, **kwargs)
            
            if result == None:
                self.ol = self.il
            
            if len(sidx) > 1: 
                del P
            if len(self.ol) == 1:
                self.ol = self.ol[0]
            if 'delete' in kwargs:
                PyRat.Data.delLayer(self.il)
            PyRat.Data.activateLayer(self.ol)
            metain = self.update(metain, PyRat.Data.shape, PyRat.Data.dtype)
            PyRat.Data.setAnnotation(metain, layer=self.ol)
            self.post()
            return self.ol
        else:
            return False
    
    def pre(self):
        pass
    
    def post(self):
        pass

    def update(self, meta, shape, dtype):
        
        if not isinstance(meta, list):
            meta  = [meta]
        ndata = len(shape) if isinstance(shape, list) else 1 
        if len(meta) != ndata:
            meta  = [meta[0]]*ndata
        if ndata == 1:
            meta = meta[0]
        return meta
        
    def read_block(self, k, idx):
        layers = self.il
        if not isinstance(layers,tuple):
            layers = (layers,)
        out = []
        for layer in layers:
            out.append(PyRat.Data.getData(block=(idx[k],idx[k]+self.blocksize,0,0),layer=layer))
        if len(layers) == 1: 
            return out[0]
        else:
            return out
    
    def save_block(self, data, k, idx, *args, **kwargs):
        if isinstance(data,tuple):
            data = list(data)
        if not isinstance(data,list):
            data = [data,]
        
        if k==0:  # first block
            self.ol = [0]*len(data)
            for d, dat in enumerate(data):
                lshape = ()
                shape  = dat.shape
                if len(shape) == 3: 
                    lshape = (shape[0],) 
                elif len(shape) == 4: 
                    lshape = (shape[0],shape[1])
                dshape = shape[-2:]
                if self.blockprocess:  dshape = (self.ydim, shape[-1])
                self.ol[d] = PyRat.Data.addLayer(dtype=dat.dtype, shape=lshape+dshape)
                if self.blockprocess:  
                    PyRat.Data.setData(dat[...,0:self.blocksize-self.blockoverlap,:],block=(idx[k],idx[k]+self.blocksize-self.blockoverlap,0,0),layer=self.ol[d])
                else:
                    PyRat.Data.setData(dat,layer=self.ol[d])
        
        elif k==len(idx)-1: # last block
            for d, dat in enumerate(data): 
                PyRat.Data.setData(dat[...,self.blockoverlap:self.blocksize,:],block=(idx[k]+self.blockoverlap,idx[k]+self.blocksize,0,0),layer=self.ol[d])
        else:    # middle block
            for d, dat in enumerate(data): 
                PyRat.Data.setData(dat[...,self.blockoverlap:self.blocksize-self.blockoverlap,:],block=(idx[k]+self.blockoverlap,idx[k]+self.blocksize-self.blockoverlap,0,0),layer=self.ol[d])
    
    def checkinput(self):
        if self.allowed_ndim != False:
             if PyRat.Data.ndim not in self.allowed_ndim:
                logging.error(self.name + ' layer dimensionality mismatch')
                return False
           
        if self.allowed_dtype != False:
            if len(set(self.allowed_dtype).intersection(PyRat.Data.dtype)) == 0:
                logging.error(self.name + ' data type mismatch')
                return False
            
        if self.require_para != False:
            annotation = PyRat.Data.getAnnotation()
            if not set(annotation.keys()).issuperset(self.require_para):
                logging.error(self.name + ' parameters missing')
                return False
        return True
    
    def calc_blocks(self, n):
        if self.blockprocess==False:
            self.blocksize = n
            self.blockoverlap = 0
            self.nthreads = 1
            return [0]
        else:
            block_start = [0]                
            block_end   = [self.blocksize]
            while block_end[-1] < n:
                    block_start.append(block_end[-1] - 2*self.blockoverlap)
                    block_end.append(block_start[-1] + self.blocksize)
            offset = block_end[-1] - n 
            block_start[-1] -= offset        
            block_end[-1]   -= offset
            return block_start
        
    def help(self):
        print self.__doc__
