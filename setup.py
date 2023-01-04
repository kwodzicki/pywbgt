import os, shutil, importlib
from setuptools import setup, convert_path, find_packages, Extension
from setuptools.command.install import install
try:
  from Cython.Build import cythonize, build_ext
except:
  raise Exception("Cython is not installed. Try installing using 'pip3 install cython'")

try:
  import numpy
except:
  raise Exception("numpy is not installed. Try installing using 'pip3 install numpy'")

NAME = 'seus_hvi_wbgt'
DESC = 'Package for SouthEast US Heat Vulnerability Index using Wet Bulb Globe Temperature.'

main_ns  = {}
ver_path = convert_path("{}/version.py".format(NAME))
with open(ver_path) as ver_file:
  exec(ver_file.read(), main_ns);

exts = [
  Extension( 
    '{}.wbgt.liljegren'.format(NAME),
    sources            = [
      '{}/wbgt/liljegren.pyx'.format(NAME),
      '{}/wbgt/src/spa.c'.format(NAME),
    ],
    extra_compile_args = ['-fopenmp'],
    #extra_link_args    = ['-lomp'],
    extra_link_args    = ['-fopenmp'],
    include_dirs       = ['{}/wbgt/src'.format(NAME)],
    depends            = ['spa.h']
  ),
  Extension( '{}.wbgt.spa'.format(NAME),
    sources            = ['{}/wbgt/spa.pyx'.format(NAME)],
    extra_compile_args = ['-fopenmp'],
    #extra_link_args    = ['-lomp'],
    extra_link_args    = ['-fopenmp'],
    include_dirs       = ['{}/wbgt/src'.format(NAME)]
  ),
  Extension( '{}.wbgt.iribarne'.format(NAME),
    sources            = ['{}/wbgt/iribarne.pyx'.format(NAME)],
    extra_compile_args = ['-fopenmp'],
    #extra_link_args    = ['-lomp'],
    extra_link_args    = ['-fopenmp'],
  )
]

setup(
  name                 = NAME,
  description          = DESC,
  url                  = "",
  author               = "",
  author_email         = "",
  version              = main_ns['__version__'],
  packages             = find_packages(),
  install_requires     = [ 'metpy' ],
  ext_modules          = cythonize( exts ),
  include_dirs         = [numpy.get_include()],
  #cmdclass             = {'build_ext' : build_ext},
  scripts              = ['bin/wbgt_gui'],
  zip_safe             = False,
)
