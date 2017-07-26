#!/usr/bin/env python
from distutils.core import setup
from Cython.Build import cythonize
import numpy
import os
import pyrat.lib.nlsar.nlsetup as nlsetup

os.environ['CFLAGS'] = '-Wno-maybe-uninitialized -Wno-cpp'

# extract internal version number (hehe, what a hack!)
pyrat_version = 'unknown'
with open('pyrat/__init__.py') as f:
    for line in f:
        if line.startswith('__version__'):
            _, _, pyrat_version = line.replace("'", '').split()
            break

# the actual setup process
setup(name='PyRAT',
      version=pyrat_version,
      description='PyRAT - Radar Tools',
      ext_modules=cythonize(["pyrat/lib/ste/*.pyx", "pyrat/filter/*.pyx", "pyrat/viewer/*.pyx"]) + [nlsetup.nlsar],
      include_dirs=[numpy.get_include()],
      packages=['pyrat', 'pyrat.save', 'pyrat.filter', 'pyrat.load', 'pyrat.filter', 'pyrat.viewer',
                'pyrat.layer', 'pyrat.insar', 'pyrat.transform', 'pyrat.polar', 'pyrat.plugins',
                'pyrat.lib', 'pyrat.lib.ste', 'pyrat.lib.templates'],
      package_data={'pyrat.lib': ['templates/*.tpl']},
      scripts=['pyrat.run'],
      data_files=[('.', ['README.md', 'requirements.txt', 'requirements_win.txt', 'MANIFEST.in', 'LICENSE',
                         'ANNOTATIONS.txt'])]
      )
