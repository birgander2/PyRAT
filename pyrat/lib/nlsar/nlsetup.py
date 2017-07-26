from distutils.core import setup, Extension
import os
import numpy as np

_lib_dir = os.path.dirname(os.path.abspath(__file__))
_lib_src = ['tools/sarerror.c', 'tools/sarprintf.c', 'tools/sarwaitbar.c', 'tools/matrixtools.c', 'tools/inv_func.c',
              'data/iorat.c', 'data/iobin.c', 'data/ioxima.c', 'data/iosar.c', 'data/ionetpbm.c',
              'data/sardata.c', 'data/rgbdata.c', 'data/fltdata.c', 'data/sar2rgb.c', 'algos/carfilter/sarboxcar.c',
              'algos/carfilter/sardiskcar.c', 'algos/carfilter/sargausscar.c', 'algos/carfilter/fltboxcar.c',
              'algos/nlsar/nlsar.c', 'algos/nlsar/sarsimstats.c', 'algos/nlsar/sarsim_glrwishart.c',
              'algos/nlsar/phi.c', 'algos/noisegen/noisegen.c']

_py_dir = os.path.join(_lib_dir,'interfaces','python')
_py_src = ['nlsartoolbox.c','sarprintf.c','sarwaitbar.c']

_src_files = [os.path.join(_lib_dir,s) for s in _lib_src] + [os.path.join(_py_dir, s) for s in _py_src]

#-O3 -W -Wall -g -fexceptions -fPIC -I/usr/include -DLAPACK -I/usr/include -I/home/jaeg_mc/git/PyRAT/pyrat/lib/nlsar/. -DOMP -fopenmp
nlsar = Extension('pyrat.lib.nlsar.nlsartoolbox',
                  language='c++',
                  define_macros = [('FFTW3F', '1'),('LAPACK', '1'),('OMP', '1')],
                  include_dirs = ['/usr/include', _lib_dir, np.get_include()],
                  extra_compile_args=['-O3', '-W', '-Wall', '-g', '-fexceptions', '-fPIC', '-fopenmp', '-Wno-error=declaration-after-statement'],
                  libraries=['m','pthread','fftw3f','lapack'],
                  library_dirs=['/usr/lib64'],
                  extra_link_args=['-fopenmp'],
                  sources =_src_files)

#setup (name = 'nlsar_test',
#       version = '1.0',
#       description = 'This is a demo package',
#       author = 'Martin v. Loewis',
#       author_email = 'martin@v.loewis.de',
#       url = 'https://docs.python.org/extending/building',
#       ext_modules = [nlsar])
