# PyRat __init__
from .cli import *
from .Worker import *
from .LayerData import *
from .FilterWorker import *
from .ImportWorker import *
from .ExportWorker import *
from .WizardWorker import *
from . import filter
from . import load
from . import save
from . import transform
from . import polar
from . import insar
from . import viewer

import os, logging, atexit, tempfile, sys
import multiprocessing

data = False
pool = False
version = 0.2


def pyrat_init(tmpdir=False, debug=False, nthreads=min(multiprocessing.cpu_count(), 8)):
    global data, pool
    if debug is True:
        logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='  %(message)s', level=logging.INFO)

    if sys.version < "3":
        logging.warning("You are running Python "+sys.version[0:3]+": Python 3.x is recommended to run PyRat!!!")
        logging.warning("Under Python2, funny things might happen. You have been warned!")

    logging.info('\n  Welcome to PyRAT (v%s)' % (version))
    logging.info('OS detected : ' + sys.platform)
    if tmpdir is False:
        if sys.platform.startswith('win'):
            configfile = os.path.join(os.path.expanduser('~'), 'pyrat.ini')
        else:
            configfile = os.path.join(os.path.expanduser('~'), '.pyratrc')
        if os.path.isfile(configfile):
            logging.debug('Found config file : ' + configfile)
            lun = open(configfile, 'rb')
            tmpdir = lun.read().rstrip()
            lun.close()
        else:
            tmpdir = tempfile.gettempdir()
            logging.info('No config file found!')
    tmpdir = tmpdir.decode()
    logging.info("Temporary directory: " + str(tmpdir))
    data = LayerData(tmpdir)
    pool = multiprocessing.Pool(nthreads)
    logging.info("Pool with " + str(nthreads) + " workers initialised")
    atexit.register(pyrat_exit)


def pyrat_exit():
    global data, pool
    logging.info('  Exiting PyRat')
    pool.close()
    for layer in data.layers.values():
        if layer.attrs['_type'] == 'Disc':
            if os.path.exists(layer.fn):
                layer.group.close()
                os.remove(layer.fn)
            logging.info('Deleting temporary file: ' + layer.fn)

interpreter = False
import inspect
inp = [s[1] for s in inspect.stack()]
try:
    if sys.ps1: interpreter = True
except AttributeError:
    if sys.flags.interactive:
        interpreter = True
if interpreter is True:
    pyrat_init()

from . import plugins
