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
    f'{NAME}.wbgt.liljegren',
    sources            = [
      f'{NAME}/wbgt/liljegren.pyx',
      f'{NAME}/wbgt/src/spa_c.c',
    ],
    extra_compile_args = [ '-fopenmp'],
    extra_link_args    = [ '-fopenmp'],
    include_dirs       = [f'{NAME}/wbgt/src'],
    depends            = [ 'spa_c.h']
  ),
  Extension( f'{NAME}.wbgt.spa',
    sources            = [f'{NAME}/wbgt/spa.pyx'],
    extra_compile_args = [ '-fopenmp'],
    extra_link_args    = [ '-fopenmp'],
    include_dirs       = [ '{}/wbgt/src'.format(NAME)]
  ),
  Extension( f'{NAME}.wbgt.iribarne',
    sources            = [f'{NAME}/wbgt/iribarne.pyx'],
    extra_compile_args = [ '-fopenmp'],
    extra_link_args    = [ '-fopenmp'],
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
