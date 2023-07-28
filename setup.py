import os
from setuptools import setup, find_packages, Extension

#from Cython.Build import cythonize
#from Cython.Build import build_ext
#cythonize
import numpy

EXTS_KWARGS = dict( 
    extra_compile_args = ['-fopenmp'],
    extra_link_args    = ['-fopenmp'],
    define_macros      = [
        ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
    ],
)

NAME = 'seus_hvi_wbgt'

# Define external extensions that must be compiled (C/Cython code)
EXTS = [
    Extension( 
        f'{NAME}.wbgt.liljegren',
        sources = [
            os.path.join('src', NAME, 'wbgt', 'liljegren.pyx'),
            os.path.join('src', NAME, 'wbgt', 'src', 'liljegren_c.c'),
            os.path.join('src', NAME, 'wbgt', 'src', 'spa_c.c'),
        ],
        include_dirs = [os.path.join('src', NAME, 'wbgt', 'src')],
        depends = ['spa_c.h'],
        **EXTS_KWARGS,
    ),
    Extension( 
        f'{NAME}.wbgt.spa',
        sources = [
            os.path.join('src', NAME, 'wbgt', 'spa.pyx'),
            os.path.join('src', NAME, 'wbgt', 'src', 'spa_c.c'),
        ],
        include_dirs = [os.path.join('src', NAME, 'wbgt', 'src')],
        depends = ['spa_c.h'],
        **EXTS_KWARGS,
    ),
    Extension( 
        f'{NAME}.wbgt.bernard',
        sources      = [os.path.join('src', NAME, 'wbgt', 'bernard.pyx')],
        **EXTS_KWARGS,
    ),
    Extension( 
        f'{NAME}.wbgt.psychrometric_wetbulb',
        sources      = [os.path.join('src', NAME, 'wbgt', 'psychrometric_wetbulb.pyx')],
        **EXTS_KWARGS,
    ),
]

# Actually run the install
setup(
    name                 = NAME,
    ext_modules          = EXTS,
#    cmdclass             = {'build_ext' : build_ext},
    include_dirs         = [numpy.get_include()],
)
