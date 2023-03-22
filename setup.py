import os
from setuptools import setup, convert_path, find_packages, Extension

try:
  from Cython.Build import cythonize
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

EXTS = [
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
  ),
]

SETUP_REQUIRES   = [
  'cython',
  'numpy',
]

INSTALL_REQUIRES = [
  'numpy',
  'pandas',
  'metpy',
  'pvlib',
]

if __name__ == "__main__":
  setup(
    name                 = NAME,
    description          = DESC,
    url                  = "",
    author               = "",
    author_email         = "",
    version              = main_ns['__version__'],
    packages             = find_packages(),
    setup_requires       = SETUP_REQUIRES,
    install_requires     = INSTALL_REQUIRES,
    include_package_data = True,
    ext_modules          = cythonize( EXTS ),
    include_dirs         = [numpy.get_include()],
    scripts              = ['bin/wbgt_gui'],
    zip_safe             = False,
  )
