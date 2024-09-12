#!/usr/bin/env python
from distutils.core import setup
from distutils.extension import Extension
import numpy
import os
import pyrat.lib.nlsar.nlsetup as nlsetup

try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

os.environ['CFLAGS'] = '-Wno-maybe-uninitialized -Wno-cpp'

# extract internal version number (hehe, what a trick!)
pyrat_version = 'unknown'
with open('pyrat/__init__.py') as f:
    for line in f:
        if line.startswith('__version__'):
            _, _, pyrat_version = line.replace("'", '').split()
            break

# set extensions (with or without cython)
ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension("pyrat.lib.ste.interpolation_extensions", ["pyrat/lib/ste/interpolation_extensions"+ext]),
              Extension("pyrat.filter.Despeckle_extensions", ["pyrat/filter/Despeckle_extensions"+ext]),
              Extension("pyrat.viewer.tools_extensions", ["pyrat/viewer/tools_extensions"+ext]),
              nlsetup.nlsar]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, force=True)
print("VERSION", pyrat_version)
# the actual setup process
setup(name='PyRAT',
      version=pyrat_version,
      description='PyRAT - Radar Tools',
      ext_modules=extensions,
      include_dirs=[numpy.get_include()],
      packages=['pyrat', 'pyrat.save', 'pyrat.filter', 'pyrat.load', 'pyrat.filter', 'pyrat.viewer',
                'pyrat.layer', 'pyrat.insar', 'pyrat.transform', 'pyrat.polar', 'pyrat.plugins',
                'pyrat.lib', 'pyrat.lib.ste', 'pyrat.lib.templates'],
      package_data={'pyrat.lib': ['templates/*.tpl']},
      scripts=['pyrat.run'],
      data_files=[('.', ['README.md', 'requirements.yml', 'MANIFEST.in'])]
      )
