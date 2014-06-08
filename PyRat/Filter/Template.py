import PyRat
import numpy as np
import pdb, logging

class Template(PyRat.FilterWorker):
    """
    Simple Template filter...
    """

    def __init__(self, *args, **kwargs):
        super(Template, self).__init__(*args, **kwargs)    
        self.name = "TEMPLATE FILTER"
        if 'win' not in self.__dict__: self.win = 7
        self.require_para  = ['PRF','CH_pol']
                
    def filter(self, array, *args, **kwargs):
        attrs = kwargs['attrs']
        track = kwargs['track']
        
        logging.error(self.name + ': I am not doing anything')
        attrs['PRF'] = 2000.0
        attrs['CH_pol'][2]  = 'RR'
        
        return array
