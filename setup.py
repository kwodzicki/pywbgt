"""
Build configuration for the pywbgt C/Cython extensions.

WHY THIS FILE EXISTS
    pywbgt ships several compiled extension modules (liljegren, bernard,
    psychrometric_wetbulb) that use OpenMP for parallelism. OpenMP is enabled
    with DIFFERENT compiler/linker flags on each platform, and the compiler is
    only known once the build actually starts. The custom ``OpenMPBuildExt``
    command below therefore injects the correct flags at build time -- this is
    what makes cross-platform binary wheels (Linux / macOS / Windows) possible.

        * Linux / GCC / Clang : -fopenmp                 (compile + link)
        * macOS / Apple clang : -Xpreprocessor -fopenmp  + link -lomp,
                                using Homebrew's libomp   (brew install libomp)
        * Windows / MSVC      : /openmp                   (no extra link flag)

    NOTE for local macOS development: historically this package was built with
    Homebrew GCC (see the archived install.sh). That still works. CI instead
    uses Apple clang + libomp, which is the standard cibuildwheel path.

CYTHON IS OPTIONAL FOR END USERS
    By DEFAULT the build compiles the pre-generated .c files committed to the
    repository, so a source build needs only a C compiler + OpenMP -- NOT
    Cython. Cython is only used to regenerate that C from the .pyx sources when
    explicitly requested via the PYWBGT_CYTHONIZE environment variable.

    Install modes (see PUBLISHING.md for detail):
        pip install pywbgt
            -> uses a prebuilt wheel if one exists (no build at all);
               otherwise builds from the committed .c files.
        pip install --no-binary pywbgt
            -> forces a source build from the committed .c files.
        PYWBGT_CYTHONIZE=1 pip install --no-binary pywbgt
            -> forces a source build that regenerates C from the .pyx sources.
               (This is also what CI sets so published wheels track the .pyx.)

    IMPORTANT: keep the committed *.c files current. After editing any .pyx,
    regenerate them with:  PYWBGT_CYTHONIZE=1 python -m build  (or
    `cythonize src/pywbgt/*.pyx`) and commit the result.
"""

import os
import sys
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import numpy

NAME = "pywbgt"

# Default to the committed .c files; regenerate from .pyx only when asked.
USE_CYTHON = os.environ.get("PYWBGT_CYTHONIZE", "0") == "1"
SRC_EXT = ".pyx" if USE_CYTHON else ".c"


def _macos_libomp_prefix():
    """Locate Homebrew's libomp (installed via ``brew install libomp``)."""
    for guess in ("/opt/homebrew/opt/libomp", "/usr/local/opt/libomp"):
        if os.path.isdir(guess):
            return guess
    try:
        return subprocess.check_output(
            ["brew", "--prefix", "libomp"], text=True
        ).strip()
    except Exception:
        return None


class OpenMPBuildExt(build_ext):
    """Inject the platform-correct OpenMP flags once the compiler is known."""

    def build_extensions(self):
        ctype = self.compiler.compiler_type
        for ext in self.extensions:
            if ctype == "msvc":
                ext.extra_compile_args += ["/openmp"]
            elif sys.platform == "darwin":
                ext.extra_compile_args += ["-Xpreprocessor", "-fopenmp"]
                ext.extra_link_args += ["-lomp"]
                prefix = _macos_libomp_prefix()
                if prefix:
                    ext.include_dirs.append(os.path.join(prefix, "include"))
                    ext.library_dirs.append(os.path.join(prefix, "lib"))
            else:  # linux / generic gcc / clang
                ext.extra_compile_args += ["-fopenmp"]
                ext.extra_link_args += ["-fopenmp"]
        super().build_extensions()


# Macros/includes shared by every extension.
DEFINE_MACROS = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
NUMPY_INCLUDE = numpy.get_include()


def _src(module):
    return os.path.join("src", NAME, module + SRC_EXT)


# liljegren also compiles against the hand-written C in src/pywbgt/src/.
EXTENSIONS = [
    Extension(
        f"{NAME}.liljegren",
        [_src("liljegren")],
        include_dirs=[os.path.join("src", NAME, "src"), NUMPY_INCLUDE],
        define_macros=DEFINE_MACROS,
    ),
    Extension(
        f"{NAME}.bernard",
        [_src("bernard")],
        include_dirs=[NUMPY_INCLUDE],
        define_macros=DEFINE_MACROS,
    ),
    Extension(
        f"{NAME}.psychrometric_wetbulb",
        [_src("psychrometric_wetbulb")],
        include_dirs=[NUMPY_INCLUDE],
        define_macros=DEFINE_MACROS,
    ),
    # Extension(f"{NAME}.solar", [_src("solar")], ...),  # disabled: solar.py is used
]

if USE_CYTHON:
    from Cython.Build import cythonize
    EXTENSIONS = cythonize(EXTENSIONS, language_level="3")

setup(
    name=NAME,
    ext_modules=EXTENSIONS,
    cmdclass={"build_ext": OpenMPBuildExt},
)
