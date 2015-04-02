#!/usr/bin/env python
from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(name='PyRAT',
      version='0.2',
      description='PyRAT - Radar Tools',
      ext_modules=cythonize(["pyrat/filter/*.pyx", "pyrat/viewer/*.pyx"]),
      include_dirs=[numpy.get_include()],
      packages=['pyrat', 'pyrat.save', 'pyrat.filter', 'pyrat.load', 'pyrat.filter', 'pyrat.tomo', 'pyrat.viewer',
                'pyrat.layer', 'pyrat.insar', 'pyrat.transform', 'pyrat.polar', 'pyrat.plugins', 'plugins'],
      scripts=['pyrat.py'],
      data_files=[('.', ['README.txt', 'requirements.txt', 'MANIFEST.in'])])
