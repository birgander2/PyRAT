# Load __init__
from __future__ import print_function

from .RatFormat import *
from .ESAR import *
from .ESAR_track import *
from .FSAR_slc import *
from .KOMPSAT5 import *
from .FSAR_dem import *
from .FSAR_track import *
from .Radarsat2 import *
from .RolfFormat import *
from .Gdal import *
from .TSX import *
from .Pixmap import *


def info():
    import sys
    from inspect import getmembers, isclass
    current_module = sys.modules[__name__]
    modules = getmembers(current_module, isclass)
    print()
    print("Content of module "+__name__+":")
    print()
    for mod in modules:
        if 'pyrat' in mod[1].__module__:
            doc = str(mod[1].__doc__)
            if doc != 'None':
                doc = doc.split('\n')[1]
            print(mod[0].ljust(20)+doc)

