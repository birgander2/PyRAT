import PyRat
import logging, pdb

class ExportWorker(object):
    def __init__(self, *args, **kwargs):
        for (k, v) in kwargs.items():         # copy keywords to self
            setattr(self, k, v)
        self.name = 'UNKNOWN'
        
    def run(self, *args, **kwargs):
        logging.info(self.name + '  '+str(dict((k, v) for k,v in self.__dict__.items() if k not in ['name','blockprocessing'])))
        self.writer(PyRat.Data.getData(), *args, **kwargs)
        return None
        
    def writer(self, array, *args, **kwargs):
        return False

