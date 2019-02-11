__version__ = '0.60-oss'

import logging, sys, os

os.environ['OMP_NUM_THREADS'] = '1'

sysexcepthook = sys.excepthook
from PIL import Image  # workaround for problems with pillow / gdal not working nicely together

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

import matplotlib.pyplot      # pseudo-import to suppress some unnecessary debug code
logging.basicConfig(format='  %(levelname)s: %(message)s', level=logging.DEBUG)

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
import json

data = False
_debug = False
_cfg = {}


# def pyrat_init(tmpdir=None, debug=False, nthreads=min(multiprocessing.cpu_count(), 4)):
def pyrat_init(debug=False, **kwargs):
    global data

    global _debug, _nthreads, _cfg
    _debug = debug

    # read config file (~/.pyratrc or the win version)

    _cfg = read_config_file(verbose=debug)

    # set number of threads

    if 'nthreads' in kwargs:
        _nthreads = kwargs['nthreads']
    elif _cfg['nthreads'] != '':
        _nthreads = int(_cfg['nthreads'])
    else:
        _nthreads = multiprocessing.cpu_count()

    # config debug mode

    pyrat_debug(debug)

    # import plugins

    import_plugins(plugin_paths=_cfg["plugindirs"], verbose=debug)
    pyrat.plugins.__name__ = "pyrat.plugins"
    pyrat.plugins.__module__ = "pyrat.plugins"
    pyrat.plugins.help = pyrat.pyrat_help("plugins", "\n  Various PyRat plugins (this can be anything!)")


    logging.info('\n  Welcome to PyRAT (v%s)' % (__version__))
    logging.info('OS detected : ' + sys.platform)

    # set up tmp dir

    if 'tmpdir' in kwargs:
        tmpdir = kwargs['tmpdir']
    elif _cfg["tmpdir"] != '':
        tmpdir = _cfg["tmpdir"]
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
    logging.info("Pool with " + str(_nthreads) + " workers initialised" + '\n')
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
        cfg["nthreads"] = cfgp["pyrat"]["nthreads"]
        cfg["plugindirs"] = [dir.strip() for dir in cfgp["pyrat"]["plugindirs"].split(',')]

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
                    logging.info(bcolors.FAIL + "Unable to import the code in plugin: %s" % filename + bcolors.ENDC)
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

