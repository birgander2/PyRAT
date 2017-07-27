__version__ = '0.4.5-oss'

import logging, sys
logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)
import scipy.misc                           # workaround for problems with pillow / gdal not working nicely together

if sys.version < "3":
    logging.error("You are running Python " + sys.version[0:3] + ": Python 3.x is required to run PyRAT!!!")
    sys.exit()
try:
    import PyQt5
except ImportError:
    logging.error("PyQt5 is required to run PyRAT. Please install requirements (see README.md) and try again !")
    sys.exit()

try:
    from osgeo import gdal
    gdal.UseExceptions()
except ImportError:
    logging.error("GDAL is required to run PyRAT. Please install requirements (see README.md) and try again !")
    sys.exit()

# # extract svn version number
# try:
#     import subprocess
#     subprocess.check_output(['svnversion', '-n']).decode()
# except:
#     __svnversion__ = 'unknown'


def docstringfrom(fromclass):
    def _decorator(func):
        func.__doc__ = fromclass.__doc__
        return func
    return _decorator


def pyrat_help(modulename, docstring):
    def f():
        import sys
        from inspect import getmembers, isclass, isfunction
        if 'plugins' in modulename:
            current_module = pyrat.plugins
            modules = getmembers(current_module, isclass)
            for mod in modules:
                mod[1].__module__ = 'pyrat.plugins.' + mod[1].__module__
        else:
            current_module = sys.modules[modulename]
            modules = getmembers(current_module, isclass)
        logging.info("")
        logging.info("Functions within the module " + modulename + ":")
        logging.info("")
        for mod in modules:
            if 'pyrat' in mod[1].__module__ and hasattr(current_module, mod[1].__name__.lower()):
                doc = str(getattr(current_module, mod[1].__name__.lower()).__doc__)
                if doc != 'None':
                    doc = doc.split('\n')[1].lstrip()
                else:
                    doc = "-"
                logging.info((mod[0].lower() + "()").ljust(30) + doc)
        logging.info("")
        logging.info(
            "Use help(" + modulename.split(".", 1)[-1] + ".function_name) for further documentation "
                                                         "on individual functions")
        logging.info("")

    f.__doc__ = docstring
    return f


from .cli import *
from .Worker import *
from .LayerData import *
from .FilterWorker import *
from .ImportWorker import *
from .ExportWorker import *
from .GroupWorker import *
from .LayerWorker import *
from . import layer
from . import filter
from . import load
from . import save
from . import transform
from . import polar
from . import insar
from . import viewer

from .tools import bcolors
import os, logging, atexit, tempfile, sys
import multiprocessing
from configparser import ConfigParser
import json

data = False
pool = False
_debug = False


def pyrat_init(tmpdir=None, debug=False, nthreads=min(multiprocessing.cpu_count(), 8)):
    global data, pool

    global _debug
    _debug = debug

    # read config file (~/.pyratrc or the win version)
    cfg = read_config_file(verbose=debug)

    # import plugins
    import_plugins(plugin_paths=cfg["plugin_paths"], verbose=debug)
    pyrat.plugins.__name__ = "pyrat.plugins"
    pyrat.plugins.__module__ = "pyrat.plugins"
    pyrat.plugins.help = pyrat.pyrat_help("plugins", "\n  Various PyRat plugins (this can be anything!)")

    pool = multiprocessing.Pool(nthreads)
    # if sys.platform.startswith('win'):
    #     for res in pool.imap(foo, [None] * nthreads):  # Workaround for delayed worker initialisation on Windows
    #         pass                                       # COMMENT: seems not to be necessary anymore !?

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    if debug is True:
        logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)
        from PyQt5 import QtCore
        QtCore.pyqtRemoveInputHook()
    else:
        logging.basicConfig(format='  %(message)s', level=logging.INFO)

    logging.info('\n  Welcome to PyRAT (v%s)' % (__version__))
    logging.info('OS detected : ' + sys.platform)

    # set up tmp dir
    if not tmpdir:
        if cfg["tmpdir"] is not None:
            tmpdir = cfg["tmpdir"]
        else:
            tmpdir = tempfile.gettempdir()
            logging.warning(
                bcolors.FAIL + bcolors.BOLD + "WARNING: Temporary directory not configured, using system default.")
            logging.warning("This often causes problems, better set it in ~/.pyratrc" + bcolors.ENDC)
    if not os.path.exists(tmpdir):
        if os.path.exists(os.path.dirname(tmpdir)):
            os.mkdir(tmpdir)
        else:
            logging.warning("WARNING: Temporary directory doesn't exist: " + tmpdir)
    logging.info("Temporary directory: " + str(tmpdir))

    data = LayerData(tmpdir)
    # pool = multiprocessing.Pool(nthreads)
    logging.info("Pool with " + str(nthreads) + " workers initialised" + '\n')
    logging.info("help() will show a list of available commands!")
    logging.info("")
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
    cfg = {}
    if config_file is None:
        if sys.platform.startswith('win'):
            config_file = os.path.join(os.path.expanduser('~'), 'pyrat.ini')
        else:
            config_file = os.path.join(os.path.expanduser('~'), '.pyratrc')

    if not os.path.isfile(config_file):
        if verbose:
            logging.info('No config file found!')
    else:
        if verbose:
            logging.debug('Found config file : ' + config_file)

        if config_type == 'json':  # load json config file
            with open(config_file) as fid:
                try:
                    cfg = json.load(fid)
                except ValueError:  # probably old-time plain ascii file
                    cfg = read_config_file(config_file, verbose=False,
                                           config_type='ini')
        if config_type == 'ini':  # load .ini config file
            try:
                cfgp = ConfigParser()
                cfgp.read(config_file)
                cfgp = {k.lower(): v for k, v in cfgp.items()}  # make case-insensitive
                cfg = dict(cfgp["pyrat"])
                if "plugin_paths" in cfgp:
                    cfg["plugin_paths"] = [v for k, v in cfgp["plugin_paths"].items()]
            except:  # probably old-time plain ascii file
                cfg = read_config_file(config_file, verbose=False,
                                       config_type='plain')
        elif config_type == 'plain':  # load initial single line config file
            lun = open(config_file, 'rb')
            tmpdir = lun.read().rstrip().decode()
            lun.close()
            cfg = {"tmpdir": tmpdir}

    # defaults, if not provided
    if "tmpdir" not in cfg:
        cfg["tmpdir"] = None
    if "plugin_paths" not in cfg:
        cfg["plugin_paths"] = []
    return cfg


class Plugins:
    """Will contain imported plugins, similar to previous module plugins"""
    pass


plugins = Plugins()


def import_plugins(plugin_paths=[], verbose=False):
    import pyrat
    pyrat_path = pyrat.__path__[0]
    default_plugin_path = os.path.dirname(pyrat_path) + "/pyrat/plugins"

    imported = []  # to import only the first occurence
    for directory in plugin_paths + [default_plugin_path]:
        if verbose:
            logging.info("Scanning for plugins: {}".format(directory))

        if not os.path.isdir(directory):
            if verbose:
                logging.debug(" - Skip. Not a directory: {}".format(directory))
            continue

        sys.path.append(directory)
        for item in os.walk(directory):
            dirpath = item[0]
            for filename in item[2]:
                if not filename.endswith(".py") or filename == "__init__.py":
                    continue
                candidate = os.path.join(dirpath, filename)

                if filename in imported:
                    continue  # don't import if another version imported
                else:
                    imported.append(filename)

                try:
                    mod = __import__(filename.split('.py')[0],
                                     globals=globals(),
                                     locals=locals(), fromlist=['*'])
                    try:
                        attrlist = mod.__all__
                    except AttributeError:
                        attrlist = dir(mod)
                    for attr in attrlist:
                        setattr(plugins, attr, getattr(mod, attr))
                    if verbose:
                        logging.info(" + Imported external plugin: " + bcolors.OKGREEN + filename + bcolors.ENDC)
                except Exception as exvar:
                    logging.info(bcolors.FAIL + "Unable to import the code in plugin: %s" % filename + bcolors.ENDC )
                    if verbose:
                        logging.debug(exvar)


def foo(bar):
    pass


def pyrat_reset():
    global data
    logging.info('Deleting PyRat Data:')
    for layer in data.layers.values():
        if layer.attrs['_type'] == 'Disc':
            if os.path.exists(layer.fn):
                layer.group.close()
                os.remove(layer.fn)
            logging.info('   ' + layer.fn)
    data = LayerData(data.tmpdir)


def pyrat_exit():
    global data, pool
    logging.info('  Exiting PyRat')
    pool.close()
    for layer in data.layers.values():
        if layer.attrs['_type'] == 'Disc':
            if os.path.exists(layer.fn):
                layer.file.close()
                os.remove(layer.fn)
            logging.debug('Deleting temporary file: ' + layer.fn)


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
