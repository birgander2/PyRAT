# Plugin __init__
from __future__ import print_function
import os, sys

for directory in ['plugins']:   # <---- here the absolute path of a plugin-dir should be!

    if not os.path.isdir(directory):
        print("skips %s (not a directory)" % directory)
        continue
    sys.path.append(directory)
    for item in os.walk(directory):
        dirpath = item[0]
        for filename in item[2]:
            if not filename.endswith(".py"):
                continue
            candidate = os.path.join(dirpath, filename)

            try:
                mod = __import__(filename.split('.py')[0], globals=globals(), locals=locals(), fromlist=['*'])
                try:
                    attrlist = mod.__all__
                except AttributeError:
                    attrlist = dir(mod)
                for attr in attrlist:
                    globals()[attr] = getattr(mod, attr)
                # print("Imported external plugin: %s" % filename)
            except Exception:
                print("Unable to import the code in plugin: %s" % filename)
                continue
