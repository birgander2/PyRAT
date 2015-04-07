# Filter __init__
import logging

from .Template import *
from .Wizard import *

from .Despeckle import *
from .Edgedetect import  *
from .Texture import  *
from .Spectrum import  *
from .Unweight import *
from .SVA import *

from .Coregister import *
from .Transforms import *
from .Presum import *
from .Console import *

def info():
    import sys
    from inspect import getmembers, isclass
    current_module = sys.modules[__name__]
    modules = getmembers(current_module, isclass)
    logging.info("")
    logging.info("Content of module "+__name__+":")
    logging.info("")
    for mod in modules:
        if 'pyrat' in mod[1].__module__:
            doc = str(mod[1].__doc__)
            if doc != 'None':
                doc = doc.split('\n')[1]
            logging.info(mod[0].ljust(20)+doc)

