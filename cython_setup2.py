from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Reverse complement C test2',
  ext_modules = cythonize("revcomp_c2.pyx"))
