import sys
import os

from setuptools import setup, find_packages, Extension
import numpy

EXT  = '.pyx' if 'build_ext' in sys.argv else '.c'
NAME = 'pywbgt'

EXTS_KWARGS = dict( 
    extra_compile_args = ['-fopenmp'],
    extra_link_args    = ['-fopenmp'],
    define_macros      = [
        ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
    ],
)

# Define external extensions that must be compiled (C/Cython code)
EXT_LILJEGREN = Extension( 
    f'{NAME}.liljegren',
    [
        os.path.join('src', NAME, 'liljegren'+EXT),
        os.path.join('src', NAME, 'src', 'spa_c.c'),
    ],
    include_dirs = [os.path.join('src', NAME, 'src')],
    depends = ['spa_c.h'],
    **EXTS_KWARGS,
)

EXT_BERNARD = Extension( 
    f'{NAME}.bernard',
    sources = [os.path.join('src', NAME, 'bernard'+EXT)],
    **EXTS_KWARGS,
)

EXT_SPA = Extension( 
    f'{NAME}.spa',
    sources = [
        os.path.join('src', NAME, 'spa'+EXT),
        os.path.join('src', NAME, 'src', 'spa_c.c'),
    ],
    include_dirs = [os.path.join('src', NAME, 'src')],
    depends = ['spa_c.h'],
    **EXTS_KWARGS,
)

EXT_PSY_WETBULB = Extension( 
    f'{NAME}.psychrometric_wetbulb',
    sources = [os.path.join('src', NAME, 'psychrometric_wetbulb'+EXT)],
    **EXTS_KWARGS,
)

EXTENSIONS = [
    EXT_LILJEGREN,
    EXT_BERNARD,
    EXT_SPA,
    EXT_PSY_WETBULB,
]

if 'build_ext' in sys.argv:
    from Cython.Build import cythonize
    _ = cythonize(
        EXTENSIONS,
        language_level = "3",
    )
    sys.exit()

# Actually run the install
setup(
    name         = NAME,
    ext_modules  = EXTENSIONS,
    include_dirs = [numpy.get_include()],
)
