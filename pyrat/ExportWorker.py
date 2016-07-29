import logging

import pyrat


class ExportWorker(pyrat.Worker):
    def __init__(self, *args, **kwargs):
        super(ExportWorker, self).__init__(*args, **kwargs)
        self.nthreads = 1

    def run(self, *args, **kwargs):

        para = [foo['var'] for foo in self.para]
        self.checkpara(kwargs, para)
        logging.info(
            self.name + '  ' + str(dict((k, v) for k, v in self.__dict__.items() if k in para or k in kwargs)))

        foo = self.open(*args, **kwargs)
        if foo is None:                        # open not overloaded -> use full image import
            if 'layer' in kwargs:
                self.writer(pyrat.data.getData(layer=kwargs['layer']), *args, **kwargs)
            else:
                self.writer(pyrat.data.getData(), *args, **kwargs)
            return self.layer
        elif foo is False:                     # open failed in some sense -> return False
            return False
        else:                                  # blockwise export
            self.layer_extract(self.block_writer, silent=False, **kwargs)
            self.close(*args, **kwargs)
            return self.layer

    def writer(self, array, *args, **kwargs):
        return False

    def block_writer(self, array, *args, **kwargs):
        return None

    def open(self, *args, **kwargs):
        return None

    def close(self, *args, **kwargs):
        return None

