import PyRat
import STEtools as STE
import numpy as np

import ipdb; stop = ipdb.set_trace

class Pixmap(PyRat.ExportWorker):
    def __init__(self, *args, **kwargs):
        super(Pixmap, self).__init__(*args, **kwargs)    
        self.name = "EXPORT TO PIXMAP"
        if 'palette' not in self.__dict__: self.palette = 0
        if 'scale' not in self.__dict__: self.scale   = 'SAR'
        if 'chscl' not in self.__dict__: self.chscl   = True
        if 'order' not in self.__dict__: self.order   = [0,1,2]
      
    def writer(self, array, *args, **kwargs):
        STE.loadct(self.palette)
        if self.scale == 'SAR' and np.iscomplexobj(array):
            array = np.abs(array)
            
        array = np.squeeze(array)
        out = np.zeros_like(array, dtype='uint8')
        if array.ndim == 3 and self.chscl == True:
            for k in range(array.shape[0]):
                array[k,...] = array[k,...] / np.mean(array[k,...])
        
        if self.scale == 'SAR':
            out = STE.sarscale(array)
        elif self.scale == 'phase':
            out = STE.phascale(array)
        elif self.scale == 'coherence':
            out = STE.cohscale(array)
        else:
            logging.error("Scaling method unknown")
            return False
        
        if array.ndim == 3: 
            out = out[self.order,...]
        STE.write_jpg(self.filename, out)
        return True
