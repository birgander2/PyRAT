# PyRat __init__

from LayerData import *
from FilterWorker import *
from ImportWorker import *
from ExportWorker import *
from LayerWorker import *
import Layer
import Filter
import Import
import Export
import Transform
import Viewer
import InSAR

import os, logging, atexit, tempfile, sys, signal
import multiprocessing
global Data, MP_Pool

def init(tmpdir = False, debug=False, nthreads=min(multiprocessing.cpu_count(), 8)):
    global Data, MP_Pool
    if debug == True:
        logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)       
    else:
        logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.INFO)
    logging.info('Welcome to PyRAT V0.1')
    if tmpdir == False:
        configfile = os.path.expanduser('~')+'/.pyratrc'
        if os.path.isfile(configfile):
            lun = open(configfile,'rb')
            tmpdir = lun.read().rstrip()
            lun.close
        else:
            tmpdir = "/tmp"
    Data = LayerData(tmpdir)
    MP_Pool = multiprocessing.Pool(nthreads)
    logging.info("Pool with "+str(nthreads)+" workers initialised")
    logging.info("Temporary directory: "+tmpdir)

@atexit.register
def __exit(*args):
    global Data, MP_Pool
    MP_Pool.close()
    for layer in Data.layers.values():
        if layer.attrs['_type'] == 'Disc':
            if os.path.exists(layer.fn):
                layer.group.close()
                os.remove(layer.fn)
            logging.info('Deleting temporary file: '+layer.fn)

