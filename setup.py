import os
from setuptools import setup, convert_path, find_packages, Extension

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
    *SETUP_REQUIRES,
    'scipy',
    'pandas',
    'pyarrow',
    'metpy',
    'pvlib',
    'matplotlib',
]

# Get version number from version file
main_ns  = {}
ver_path = convert_path( os.path.join(NAME, 'version.py') )
with open(ver_path) as ver_file:
  exec(ver_file.read(), main_ns)

# Set up some keyword arguments for compiling external extensions
# (C/Cython code)
EXTS_KWARGS = dict( 
    extra_compile_args = ['-fopenmp'],
    extra_link_args    = ['-fopenmp'],
)

# Define external extensions that must be compiled (C/Cython code)
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

# Actually run the install
if __name__ == "__main__":
    from Cython.Build import cythonize
    import numpy

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
