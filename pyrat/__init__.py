# PyRat __init__
import logging
logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)

from .cli import *
from .Worker import *
from .LayerData import *
from .FilterWorker import *
from .ImportWorker import *
from .ExportWorker import *
from .WizardWorker import *
from .LayerWorker import *
from . import layer
from . import filter
from . import load
from . import save
from . import transform
from . import polar
from . import insar
from . import viewer

import os, logging, atexit, tempfile, sys
import multiprocessing
from configparser import ConfigParser
import json

data = False
pool = False
_debug = False
version = 0.3


def pyrat_init(tmpdir=None, debug=False, nthreads=min(multiprocessing.cpu_count(), 8)):
    global data, pool

    global _debug
    _debug = debug

    pool = multiprocessing.Pool(nthreads)
    if sys.platform.startswith('win'):
        for res in pool.imap(foo, [None]*nthreads):      # Workaround for delayed worker initialisation on Windows
            pass

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    if debug is True:
        logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='  %(message)s', level=logging.INFO)

    if sys.version < "3":
        logging.warning("You are running Python "+sys.version[0:3]+": Python 3.x is recommended to run PyRat!!!")
        logging.warning("Under Python2, funny things might happen. You have been warned!")

    logging.info('\n  Welcome to PyRAT (v%s)' % (version))
    logging.info('OS detected : ' + sys.platform)

    # read config file (~/.pyratrc or the win version)
    cfg = read_config_file()

    # set up tmp dir
    if tmpdir is None:
        if "tmpdir" in cfg:
            tmpdir = cfg["tmpdir"]
        else:
            tmpdir = tempfile.gettempdir()
    if not os.path.exists(tmpdir):
        if os.path.exists(os.path.dirname(tmpdir)):
            os.mkdir(tmpdir)
        else:
            logging.warning("WARNING: Temporary directory doesn't exist: " + tmpdir)


    logging.info("Temporary directory: " + str(tmpdir))
    data = LayerData(tmpdir)
    # pool = multiprocessing.Pool(nthreads)
    logging.info("Pool with " + str(nthreads) + " workers initialised" + '\n')
    atexit.register(pyrat_exit)


def read_config_file(config_file=None, verbose=True, config_type='json'):
    """Read config file into a dict structure.

    Three config file formats supported:
    - json (default)
    - ini files
    - plain (tmpdir only)

    If config_file is found, it tries to read it first in json format. If it
    doesn't work, it tries to read in the .ini format. If this doesn't work
    as well, it will finally attempt to read it in the plain form.

    config_type : {'json', 'ini', 'plain'}
    """
    if config_file is None:
        if sys.platform.startswith('win'):
            config_file = os.path.join(os.path.expanduser('~'), 'pyrat.ini')
        else:
            config_file = os.path.join(os.path.expanduser('~'), '.pyratrc')

    if not os.path.isfile(config_file):
        if verbose:
            logging.info('No config file found!')
        return {}

    if verbose:
        logging.debug('Found config file : ' + config_file)

    if config_type == 'json':   # load json config file
        with open(config_file) as fid:
            try:
                cfg = json.load(fid)
            except ValueError:  # probably old-time plain ascii file
                cfg = read_config_file(config_file, verbose=False,
                                       config_type='ini')
    if config_type == 'ini':   # load .ini config file
        try:
            cfgp = ConfigParser()
            cfgp.read(config_file)
            cfgp = {k.lower():v for k,v in cfgp.items()} # make case-insensitive
            cfg = cfgp["pyrat"]
        except:  # probably old-time plain ascii file
            cfg = read_config_file(config_file, verbose=False,
                                   config_type='plain')
    elif config_type == 'plain': # load initial single line config file
        lun = open(config_file, 'rb')
        tmpdir = lun.read().rstrip().decode()
        lun.close()
        cfg = {"tmpdir": tmpdir}

    return cfg


def foo(bar):
    pass


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
