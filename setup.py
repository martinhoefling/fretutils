from distutils.core import setup, Extension, Command
from Cython.Distutils import build_ext
import sys, os

MKLPATH = "/usr/local/intel/mkl/10.2.0.013/lib/em64t"

fretnumpyext = Extension('FRETUtils.fretnumpyext',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '1')],
#                    include_dirs = ['/usr/local/sage-4.1-intel/local/include/python2.6','/usr/local/sage-4.1-intel/local/lib/python2.6/site-packages/numpy/core/include'],
#                    include_dirs = ['/usr/include/python2.6','/usr/lib/python2.6/dist-packages/numpy/core/include'],
                    include_dirs = ['src', '/home/mhoefli/local-cluster/include', '/usr/include/python2.6', '/usr/lib64/python2.6/site-packages/numpy/core/include', '/usr/include/python2.7', '/usr/lib64/python2.7/site-packages/numpy/core/include', '/usr/local/intel/mkl/10.2.0.013/include'],
                    libraries = [],
                    library_dirs = ['/cm/shared/apps/intel/Compiler/11.1/046/lib/intel64/', '/usr/local/intel/mkl/10.2.0.013/lib/em64t'],
                    sources = ['src/fretnumpyext.c', 'src/SFMT.c'],
                    extra_compile_args = ['-D_FORTIFY_SOURCE=0', '-DMEXP=607', '-msse2', '-DHAVE_SSE2', '-DUSE_MT'],
#                    extra_link_args=['-L%s %s/libmkl_solver_ilp64_sequential.a -Wl,--start-group -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread'%(MKLPATH,MKLPATH)]
                   )

cythonext = Extension("FRETUtils.PhotonGenerator", ["lib/FRETUtils/PhotonGenerator.pyx"],
                      include_dirs = ['/usr/lib64/python2.7/site-packages/numpy/core/include'])

setup (name = 'fretutils',
       version = '1.1',
       description = 'This module provides FRET MC routines.',
       author = 'Martin Hoefling',
       author_email = 'mhoefli@gwdg.de',
       url = '',
       long_description = '''
FRET MC routines.
''',
       package_dir = {'': 'lib'},
       packages = ['FRETUtils', 'FRETUtilsTests'],
       cmdclass = {"build_ext": build_ext},
       scripts = ['scripts/md2fret.py'],
       ext_modules = [fretnumpyext, cythonext])


