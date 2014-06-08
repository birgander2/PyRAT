import PyRat
import logging, copy, pdb, time, curses, sys
import numpy as np
import multiprocessing as mp

def exec_parallel(args):
    return args[0].filter(args[1], track=args[2], meta=args[3])
 
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
        self.nthreads      = 8 # mp.cpu_count() 
        self.noshow = ['noshow','nthreads','blockoverlap','blocksize', 'name','blockprocess','allowed_dtype','allowed_ndim','require_para']
    
    def filter(self, data, *args, **kwargs):
        logging.error(self.name + ': No filter method defined')
        return False
    
    def run_single(self, *args, **kwargs):
        print "SINGLE"
        logging.info(self.name + '  '+str(dict((k, v) for k,v in self.__dict__.items() if k not in self.noshow)))
        if self.checkinput():
            meta  = PyRat.Data.getAnnotation()
            track = PyRat.Data.getTrack()
            result   = self.filter(PyRat.Data.getData(), track=track, meta=meta, *args, **kwargs)
            newlayer = PyRat.Data.addLayer(result, *args, **kwargs)
            PyRat.Data.setAnnotation(meta, layer=newlayer)
            PyRat.Data.setTrack(track, layer=newlayer)
            return newlayer
        else:
            return False
    def run(self, *args, **kwargs):
        logging.info(self.name + '  '+str(dict((k, v) for k,v in self.__dict__.items() if k not in self.noshow)))
        if self.checkinput():
            if 'replace' in kwargs and self.blockprocess==True:
                logging.warning("Blockprocessing and replace keyword not compatible - ignoring it!")
                del kwargs['replace']
            self.il   = PyRat.Data.active
            meta      = PyRat.Data.getAnnotation(layers=self.il)
            track     = PyRat.Data.getTrack(layer=self.il)
            self.ydim = PyRat.Data.dshape[0]
            if isinstance(self.ydim, tuple):
                self.ydim = self.ydim[0]
            sidx      = self.calc_blocks(self.ydim)

            if len(sidx) > 1 and self.nthreads > 1: 
                idx = [sidx[i:i+self.nthreads] for i in range(0, len(sidx), self.nthreads)]
            else:
                idx  = sidx
                
            if len(sidx) > 1:
                P = ProgressBar('  '+self.name, len(sidx)-1)
                for k in range(len(idx)):
                    metain = copy.copy(meta)
                    if self.nthreads > 1:
                        bidx = idx[k]
                        inputs = []
                        for l in range(len(bidx)):
                            data = self.read_block(l, bidx)
                            inputs.append((self,data,track,meta))
                        result = PyRat.MP_Pool.imap(exec_parallel, inputs)
                        for l, data in enumerate(result):
                            pos = sidx.index(bidx[l])
                            self.save_block(data, pos, sidx, *args, **kwargs)
                            P.update(pos)
                    else:
                        data = self.read_block(k, sidx)
                        data = self.filter(data, track=track, meta=metain, *args, **kwargs)
                        self.save_block(data, k, sidx) 
                        if len(sidx) > 1: P.update(k)
            else:
                metain = copy.copy(meta)
                result  = self.filter(PyRat.Data.getData(), track=track, meta=metain, *args, **kwargs)
                self.ol = PyRat.Data.addLayer(result, *args, **kwargs)
   
            if len(sidx) > 1: 
                del P
            if len(self.ol) == 1:
                self.ol = self.ol[0]
            
            PyRat.Data.setAnnotation(metain, layers=self.ol)
            PyRat.Data.setTrack(track, layer=self.ol)
            PyRat.Data.activateLayer(self.ol)
            return self.ol
        else:
            return False
        
    def process_block(self, k):
        pass
    
    def read_block(self, k, idx):
        layers = self.il
        if not isinstance(layers,tuple):
            layers = (layers,)
        out = []
        for layer in layers:
            out.append(PyRat.Data.getData(block=[idx[k],idx[k]+self.blocksize,0,0],layer=layer))
        if len(layers) == 1: 
            return out[0]
        else:
            return tuple(out)
    
    def save_block(self, data, k, idx, *args, **kwargs):
        if not isinstance(data,tuple):
            data = (data,)
        
        if k==0:  # first block
            self.ol = [0]*len(data)
            for d, dat in enumerate(data):
                lshape = (1,)
                shape  = dat.shape
                if len(shape) == 3: 
                    lshape = (shape[0],) 
                elif len(shape) == 4: 
                    lshape = (shape[0],shape[1])
                dshape = shape[-2:]
                if self.blockprocess:  dshape = (self.ydim, shape[-1])
                self.ol[d] = PyRat.Data.addLayer(dtype=dat.dtype, shape=lshape+dshape, *args, **kwargs)
                if self.blockprocess:  
                    PyRat.Data.setData(dat[...,0:self.blocksize-self.blockoverlap,:],block=[idx[k],idx[k]+self.blocksize-self.blockoverlap,0,0],layer=self.ol[d])
                else:
                    PyRat.Data.setData(dat,layer=self.ol[d])
            self.ol = tuple(self.ol)
        
        elif k==len(idx)-1: # last block
            for d, dat in enumerate(data): 
                PyRat.Data.setData(dat[...,self.blockoverlap:self.blocksize,:],block=[idx[k]+self.blockoverlap,idx[k]+self.blocksize,0,0],layer=self.ol[d])
        else:    # middle block
            for d, dat in enumerate(data): 
                PyRat.Data.setData(dat[...,self.blockoverlap:self.blocksize-self.blockoverlap,:],block=[idx[k]+self.blockoverlap,idx[k]+self.blocksize-self.blockoverlap,0,0],layer=self.ol[d])
    
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
        
class ProgressBar():
    """
    Simple progress bar for the command line. 
   
    :author: Andreas Reigber
    """
    def __init__(self, message, max, **kwargs):
        curses.setupterm()
        terminal_width = curses.tigetnum('cols')
        if not sys.stdout.isatty(): terminal_width = 70
        self.message = message.ljust(25)
        if terminal_width < 50: self.message = self.message[:terminal_width/2]
        self.width   = terminal_width - len(self.message) - 10
        self.max     = max

    def __del__(self):
        print
        
    def update(self,val):
        """
        Updates the progress bar to the given value of progress
        """
        percent = float(val) / self.max
        hashes = '#' * int(round(percent * self.width))
        spaces = ' ' * (self.width - len(hashes))
        retline = "\r" if sys.stdout.isatty() else ""
        if sys.stdout.isatty() or val == 0:
            sys.stdout.write(retline+self.message+": [{0}] {1}%".format(hashes + spaces, int(round(percent * 100))))
            sys.stdout.flush()
           
