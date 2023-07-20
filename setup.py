import os
from setuptools import setup, find_packages, Extension

from Cython.Build import cythonize
import numpy

EXTS_KWARGS = dict( 
    extra_compile_args = ['-fopenmp'],
    extra_link_args    = ['-fopenmp'],
)

NAME = 'seus_hvi_wbgt'

# Define external extensions that must be compiled (C/Cython code)
EXTS = [
#    Extension( 
#        f'{NAME}.wbgt.liljegren',
#        sources            = [
#            os.path.join('src', NAME, 'wbgt', 'liljegren.pyx'),
#            os.path.join('src', NAME, 'wbgt', 'src', 'spa_c.c'),
#        ],
#        include_dirs = [os.path.join('src', NAME, 'wbgt', 'src')],
#        depends      = ['spa_c.h'],
#        **EXTS_KWARGS,
#    ),
    Extension( 
        f'{NAME}.wbgt.spa',
        sources      = [
            os.path.join('src', NAME, 'wbgt', 'spa.pyx'),
            os.path.join('src', NAME, 'wbgt', 'src', 'spa_c.c'),
        ],
        include_dirs = [os.path.join('src', NAME, 'wbgt', 'src')],
        depends      = ['spa_c.h'],
        **EXTS_KWARGS,
    ),
#    Extension( 
#        f'{NAME}.wbgt.bernard',
#        sources      = [os.path.join('src', NAME, 'wbgt', 'bernard.pyx')],
#        **EXTS_KWARGS,
#    ),
#    Extension( 
#        f'{NAME}.wbgt.psychrometric_wetbulb',
#        sources      = [os.path.join('src', NAME, 'wbgt', 'psychrometric_wetbulb.pyx')],
#        **EXTS_KWARGS,
#    ),
]

# Actually run the install
setup(
    name                 = NAME,
#    packages             = find_packages(where="src"),
#    package_dir          = {"" : "src"},
#    include_package_data = True,
    ext_modules          = cythonize( EXTS ),
    include_dirs         = [numpy.get_include()],
)
