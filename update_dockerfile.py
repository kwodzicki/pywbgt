"""
Script to update Dockerfile with dependancies

This script is designed to update/replace the last line in a Dockerfile
with a 'pip install' command that contains all the packages required for
running/installing the 'main' package for the container.

This script assumes that there is a 'INSTALL_REQUIRES' global variable
in setup.py that defines all the packages required for installation.
If there are special links required to find packages, those should be
defined in the 'DEPENDENCY_LINKS' global variable of setup.py.

Please ensure that the call the 'setuptools.setup()' is wrapped in a 

    if __name__ == "__main__":

call to ensure that package installation does not begin when the 
setup.py file is imported below.

"""

import os
import sys
import argparse

if __name__ == "__main__":
    import setup

    parser = argparse.ArgumentParser(
        description='Run to update Dockerfile to include "install_requires" packages in container build',
    )
    parser.add_argument(
        'dockerfile',
        type=str,
        help='Path to Docker file to update',
    )

    args = parser.parse_args()
    dockerfile = os.path.abspath( args.dockerfile )

    if not hasattr(setup, 'INSTALL_REQUIRES'):
        sys.exit(1)
    requires = ' '.join(
        [f'"{r}"' for r in setup.INSTALL_REQUIRES]
    )

    links = (
        ' '.join( [f'-f "{l}"' for l in setup.DEPENDENCY_LINKS] )
        if hasattr(setup, 'DEPENDENCY_LINKS') else
        ''
    )

    with open(dockerfile, mode="r") as fid:
        lines = fid.read().splitlines()

    lines[-1] = f"RUN pip install --upgrade --upgrade-strategy only-if-needed {requires} {links}"

    with open(dockerfile, mode="w") as fid:
        fid.write( os.linesep.join(lines) )
