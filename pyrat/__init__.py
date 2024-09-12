__version__ = '0.65+oss'

import logging, sys, os

if 'OMP_NUM_THREADS' not in os.environ:
    os.environ['OMP_NUM_THREADS'] = '1'

sysexcepthook = sys.excepthook
from PIL import Image  # workaround for problems with pillow / gdal not working nicely together

# blacklist pooch to avoid PermissionError when importing skimage (which creates directories!)
sys.modules['pooch'] = None

if sys.version < "3":
    logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)
    logging.error("You are running Python " + sys.version[0:3] + ": Python 3.x is required to run PyRAT!!!")
    sys.exit()

try:
    import PyQt5
except ImportError:
    logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)
    logging.error("PyQt5 is required to run PyRAT. Please install requirements (see README.md) and try again !")
    sys.exit()

try:
    from osgeo import gdal
    gdal.UseExceptions()
except ImportError:
    logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)
    logging.error("GDAL is required to run PyRAT. Please install requirements (see README.md) and try again !")
    sys.exit()


# # extract svn version number
# try:
#     import subprocess
#     subprocess.check_output(['svnversion', '-n']).decode()
# except:
#     __svnversion__ = 'unknown'


def exithook(var1, var2, var3):
    if _debug is True:
        sysexcepthook(var1, var2, var3)
    else:
        import traceback
        tb = sys.exc_info()[2]
        tbinfo = traceback.extract_tb(tb)[-1]
        logging.error(bcolors.FAIL + 'ERROR : ' + str(var2))
        logging.error(var1.__name__ + " in " + os.path.basename(tbinfo[0]) +
                      " (line " + str(tbinfo[1]) + ")" + bcolors.ENDC)


def docstringfrom(fromclass):
    def _decorator(func):
        func.__doc__ = fromclass.__doc__
        return func

    return _decorator


def pyrat_help(modulename, docstring):
    def f():
        import sys
        from inspect import getmembers, isclass
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

import matplotlib.pyplot      # pseudo-import to suppress some unnecessary debug code

from . import layer
from . import filter
from . import load
from . import save
from . import transform
from . import polar
from . import insar
from . import viewer

from .tools import bcolors
import logging, atexit, tempfile, sys
import multiprocessing
from configparser import ConfigParser


# initialise global variables
data = False
_debug = False
_cfg = {}
_nthreads = multiprocessing.cpu_count()


# def pyrat_init(tmpdir=None, debug=False, nthreads=min(multiprocessing.cpu_count(), 4)):
def pyrat_init(debug=False, **kwargs):
    global data

    global _debug, _nthreads, _cfg, _sc
    _debug = debug

    # read config file (~/.pyratrc or the win version)

    _cfg = read_config_file(verbose=debug)

    # set number of threads

    if 'nthreads' in kwargs:
        _nthreads = kwargs['nthreads']
    elif _cfg['nthreads'] != '':
        _nthreads = int(_cfg['nthreads'])

    # config debug mode

    pyrat_debug(debug)

    # set up silent mode
    if 'silent' in kwargs:
        silent = kwargs['silent']
    else:
        silent = False

    # import plugins

    import_plugins(plugin_paths=_cfg["plugindirs"], verbose=debug)

    if not silent:
        logging.info('\n  Welcome to PyRAT (v%s)' % (__version__))
        logging.info('OS detected : ' + sys.platform)

    # set up tmp dir

    if 'tmpdir' in kwargs:
        tmpdir = kwargs['tmpdir']
    elif _cfg["tmpdir"] != '':
        tmpdir = _cfg["tmpdir"]
    else:
        tmpdir = tempfile.gettempdir()
        if not silent:
            logging.warning(
                bcolors.FAIL + bcolors.BOLD + "WARNING: Temporary directory not configured, using system default.")
            logging.warning("This often causes problems, better set it in ~/.pyratrc" + bcolors.ENDC)

    if not os.path.exists(tmpdir):
        if os.path.exists(os.path.dirname(tmpdir)):
            os.mkdir(tmpdir)
        else:
            if not silent:
                logging.warning("WARNING: Temporary directory doesn't exist: " + tmpdir)

    if not silent:
        logging.info("Temporary directory: " + str(tmpdir))

    if _cfg['spark'] is not None:
        # initialise the Spark context when the configuration contains a spark setting such as
        # spark: yarn
        # for processing using an entire YARN cluster (or local[8] for using just 8 local cores)
        _sc = get_pyrat_spark_conext(_cfg['spark'])
        if not silent:
            logging.info("Initialised Spark context for multi-processing" + '\n')
    else:
        _sc = None
        if not silent:
            logging.info("Pool with " + str(_nthreads) + " workers initialised" + '\n')

    data = LayerData(tmpdir)
    if not silent:
        logging.info("help() will show a list of available commands!")

    atexit.register(pyrat_exit)
    sys.excepthook = exithook


def pyrat_debug(debug):
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    if debug is True:
        logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)
        from PyQt5 import QtCore
        QtCore.pyqtRemoveInputHook()
    else:
        logging.basicConfig(format='  %(message)s', level=logging.INFO)


def read_config_file(config_file=None, verbose=True):
    """Read config file into a dict structure.

    """
    cfg = {}
    if config_file is None:
        if sys.platform.startswith('win'):
            config_file = os.path.join(os.path.expanduser('~'), 'pyrat.ini')
        else:
            config_file = os.path.join(os.path.expanduser('~'), '.pyratrc')

    if os.path.isfile(config_file):
        if verbose:
            logging.debug('Found config file : ' + config_file)

        cfgp = ConfigParser()
        try:
            cfgp.read(config_file)
        except:
            logging.warning(bcolors.FAIL + "Error reading config file (old format maybe?)" + bcolors.ENDC)
            new_conigfile(config_file)
            cfgp.read(config_file)

        cfgp = {k.lower(): v for k, v in cfgp.items()}  # make case-insensitive
        cfgp['pyrat'] = {k.lower(): v for k, v in cfgp['pyrat'].items()}

        cfg["tmpdir"] = cfgp["pyrat"]["tempdir"]
        cfg["nthreads"] = cfgp["pyrat"].get("nthreads",8)
        cfg["plugindirs"] = [dir.strip() for dir in cfgp["pyrat"]["plugindirs"].split(',')]
        cfg["spark"] = cfgp["pyrat"].get("spark", None)

        if "qgis" in cfgp:
            cfg["qgispath"] = cfgp["qgis"]["path"]

    else:
        logging.warning(bcolors.FAIL + "No config file found!"  + bcolors.ENDC)
        new_conigfile(config_file)
        cfg = read_config_file(config_file, verbose=verbose)
    return cfg


def new_conigfile(config_file):
    logging.warning(bcolors.FAIL + "Generatig a new config file: '" + config_file + "'. Please review it!" + bcolors.ENDC)
    logging.warning(bcolors.FAIL + "In particular choose a better TempDir setting than the system's default" + bcolors.ENDC)
    with open(config_file, 'w') as lun:
        lun.write("[PYRAT]" + "\n")
        lun.write("TempDir: " + tempfile.gettempdir() + "\n")
        lun.write("PluginDirs: " + "\n")
        lun.write("NThreads: " + str(multiprocessing.cpu_count()) + "\n")


def import_plugins(plugin_paths=[], verbose=False):
    import importlib.abc
    import importlib.util

    class CustomImporter(importlib.abc.Loader, importlib.abc.MetaPathFinder):

        def __init__(self, verbose: bool, plugin_paths: list):
            super().__init__()
            default_plugin_path = os.path.dirname(__path__[0]) + "/pyrat/plugins"
            self.plugin_paths = plugin_paths + [default_plugin_path]
            self.verbose = verbose

        def find_spec(self, fullname, path, *args):
            if fullname == "pyrat.plugins":
                return importlib.util.spec_from_loader(fullname, self)  # Return itself

        def exec_module(self, module):
            pathbackup = sys.path

            for directory in self.plugin_paths:
                if not os.path.isdir(directory):
                    if self.verbose:
                        logging.info(" - Skip. Not a directory: {}".format(directory))
                    continue

                sys.path = [directory]
                for item in os.walk(directory):
                    for filename in item[2]:
                        if not filename.endswith(".py"):
                            continue

                        modname = filename.split('.py')[0]
                        try:
                            spec = importlib.util.find_spec(modname)
                            spec.loader.name = 'pyrat.plugins.' + modname
                            mod = spec.loader.load_module('pyrat.plugins.' + modname)
                        except Exception as exvar:
                            logging.info(
                                bcolors.FAIL + "Unable to import the code in plugin: %s" % filename + bcolors.ENDC)
                            if self.verbose:
                                logging.debug(exvar)
                            continue

                        setattr(module, modname, mod)
                        for attrname in dir(mod):
                            if attrname.startswith("__"):  # Do not import module internals
                                continue

                            setattr(module, attrname, getattr(mod, attrname))
                        if self.verbose:
                            logging.info(" + Imported external plugin: " + bcolors.OKGREEN + filename + bcolors.ENDC)

            module.help = pyrat.pyrat_help("pyrat.plugins", "\n  Various PyRat plugins (this can be anything!)")
            sys.path = pathbackup
            return module

    sys.meta_path.insert(0, CustomImporter(verbose, plugin_paths))
    from . import plugins


def get_pyrat_spark_conext(spark_master):
    """
    initialises the PyRAT spark context, optionally restoring some conda environment settings
    :param spark_master: spark master URL, such as 'local[8]' or 'yarn'
    :return: initialised spark context
    """
    try:
        # normaly, this import should work (e.g. when running pyrat in a properly configured pyrat VE)
        import pyspark as spark
    except ImportError:
        # special case: pyrat was started with the correct python interpreter, but the environment
        # variables for the conda environment have not been set.
        # This happens when running pyrat in pycharm, where the correct python interpreter is used
        # but the actual environment is never fully "activated". In this case, we restore the conda
        # environment settings explicitly!
        conda_state = os.path.join(os.path.dirname(sys.executable), '..', 'conda-meta', 'state')
        try:
            with open(conda_state, 'r') as f:
                import json
                state = json.load(f)
                os.environ.update(state['env_vars'])
            path_ext = [p for p in os.environ.get('PYTHONPATH', '').split(':') if len(p) > 0]
            sys.path.extend(path_ext)
        except IOError:
            pass
        import pyspark as spark
    # create and return the contex
    conf = spark.SparkConf(True).setAppName('PyRAT')
    conf.setMaster(spark_master)
    return spark.SparkContext(conf=conf)


def foo(bar):
    pass


def pyrat_reset():
    global data
    logging.info('Deleting PyRat Data:')
    for layer in data.layers.values():
        if layer.attrs['_type'] == 'Disc':
            if os.path.exists(layer.fn):
                layer.file.close()
                os.remove(layer.fn)
            logging.info('   ' + layer.fn)
    data = LayerData(data.tmpdir)


def pyrat_exit():
    global data
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
    if sys.ps1:
        interpreter = True
except AttributeError:
    if sys.flags.interactive:
        interpreter = True
    # if sys.__stdin__.isatty():
    #     interpreter = True

try: # Check for QGIS python console
    import qgis.utils
    interpreter = True
except:
    pass

if interpreter is True:
    pyrat_init()

# Launched from within QGIS?
try:
    if "qgispath" in _cfg:
        sys.path.append(_cfg["qgispath"])
    from .qgis import *
except ImportError:
   pass
   # logging.error("QGis is required to run the Qgis interface!")

