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

NAME  = 'seus_hvi_wbgt'
DESC  = 'Package for SouthEast US Heat Vulnerability Index using Wet Bulb Globe Temperature.'
AUTH  = 'Kyle R. Wodzicki'
EMAIL = 'kwodzicki@cicsnc.org'
URL   = ''

SETUP_REQUIRES   = [
    'cython',
    'numpy',
]

INSTALL_REQUIRES = [
    'numpy',
    'scipy',
    'pandas',
    'metpy',
    'pvlib',
    'matplotlib',
]

main_ns  = {}
ver_path = convert_path( os.path.join(NAME, 'version.py') )
with open(ver_path) as ver_file:
  exec(ver_file.read(), main_ns)

EXTS_KWARGS = dict( 
    extra_compile_args = ['-fopenmp'],
    extra_link_args    = ['-fopenmp'],
)

EXTS = [
    Extension( 
        f'{NAME}.wbgt.liljegren',
        sources            = [
            os.path.join(NAME, 'wbgt', 'liljegren.pyx'),
            os.path.join(NAME, 'wbgt', 'src', 'spa_c.c'),
        ],
        include_dirs = [os.path.join(NAME, 'wbgt', 'src')],
        depends      = ['spa_c.h'],
        **EXTS_KWARGS,
    ),
    Extension( 
        f'{NAME}.wbgt.spa',
        sources      = [
            os.path.join(NAME, 'wbgt', 'spa.pyx'),
            os.path.join(NAME, 'wbgt', 'src', 'spa_c.c'),
        ],
        include_dirs = [os.path.join(NAME, 'wbgt', 'src')],
        depends      = ['spa_c.h'],
        **EXTS_KWARGS,
    ),
    Extension( 
        f'{NAME}.wbgt.bernard',
        sources      = [os.path.join(NAME, 'wbgt', 'bernard.pyx')],
        **EXTS_KWARGS,
    ),
    Extension( 
        f'{NAME}.wbgt.psychrometric_wetbulb',
        sources      = [os.path.join(NAME, 'wbgt', 'psychrometric_wetbulb.pyx')],
        **EXTS_KWARGS,
    ),
]

if __name__ == "__main__":
    setup(
        name                 = NAME,
        description          = DESC,
        url                  = URL,
        author               = AUTH,
        author_email         = EMAIL,
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
