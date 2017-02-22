#!/usr/bin/env python
from distutils.core import setup
from Cython.Build import cythonize
import numpy

# # write svn version number (if possible)
# try:
#     import subprocess
#     __svnversion__ = subprocess.check_output(['svnversion', '-n']).decode()
# except:
#     __svnversion__ = 'unknown'
# with open('pyrat/VERSION', 'w') as f:
#     f.write(__svnversion__)

# extract internal version number (hehe, what a trick!)

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
      ext_modules=cythonize(["pyrat/lib/ste/*.pyx", "pyrat/filter/*.pyx", "pyrat/viewer/*.pyx"]),
      include_dirs=[numpy.get_include()],
      packages=['pyrat', 'pyrat.save', 'pyrat.filter', 'pyrat.load', 'pyrat.filter', 'pyrat.viewer',
                'pyrat.layer', 'pyrat.insar', 'pyrat.transform', 'pyrat.polar', 'pyrat.plugins'],
      scripts=['PyRat'],
      data_files=[('.', ['README.txt', 'requirements.txt', 'MANIFEST.in'])]
      )
