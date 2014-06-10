import PyRat, logging
import numpy
import IDL

import ipdb; stop = ipdb.set_trace

class Coherence(PyRat.FilterWorker):
    """
    Calc InSAR phase between two images
    """

    def __init__(self, *args, **kwargs):
        super(Coherence, self).__init__(*args, **kwargs)    
        self.name = "CALC InSAR COHERENCE"
        if 'win'   not in self.__dict__: self.win=[7,7]
        self.blockoverlap = self.win[0]/2+1
        self.blockprocess = True
        self.nthreads = 1
        
    def filter(self, array, *args, **kwargs):
        coh = numpy.abs(IDL.smooth(array[0]*numpy.conj(array[1]),self.win)
                        / numpy.sqrt(IDL.smooth(array[0]*numpy.conj(array[0]),self.win)
                        * IDL.smooth(array[1]*numpy.conj(array[1]),self.win)))
        return numpy.clip(numpy.nan_to_num(coh),0.0,1.0)  # get rid of numerical inaccuracies!
    
        
    