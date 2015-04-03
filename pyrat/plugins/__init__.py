# Plugin __init__
from __future__ import print_function
import os, sys


def import_plugins(plugin_paths = None, verbose=False):
    if plugin_paths is None:
        import pyrat
        pyrat_path = pyrat.__path__[0]
        plugin_paths = ["plugins",                    # local directory
                        pyrat_path + "/plugins",      # pyrat src directory
                        pyrat_path + "/../plugins"]   # pyrat root directory

    for directory in plugin_paths:   # <---- here the absolute path of a plugin-dir should be!

        if verbose:
            print("Scanning for plugins: {}".format(directory))

        if not os.path.isdir(directory):
            if verbose:
                print(" - Skip this. Not a directory: {}".format(directory))
            continue

        sys.path.append(directory)
        for item in os.walk(directory):
            dirpath = item[0]
            for filename in item[2]:
                if not filename.endswith(".py") or filename == "__init__.py":
                    continue
                candidate = os.path.join(dirpath, filename)

                try:
                    mod = __import__(filename.split('.py')[0], globals=globals(),
                                     locals=locals(), fromlist=['*'])
                    try:
                        attrlist = mod.__all__
                    except AttributeError:
                        attrlist = dir(mod)
                    for attr in attrlist:
                        globals()[attr] = getattr(mod, attr)
                    if verbose:
                        print(" + Imported external plugin: %s" % filename)
                except Exception:
                    print("Unable to import the code in plugin: %s" % filename)
                    continue


import_plugins(verbose=True)
