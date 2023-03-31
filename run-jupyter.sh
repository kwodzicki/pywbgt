#!/bin/bash

DOCKERFILE=Dockerfile

# Replace last 2 lines of Dockerfile with packages to
# install and cmd to run when container starts
#
python3 - << EOF
import os
from setup import INSTALL_REQUIRES

requires = ' '.join(
  [f'\"{r}\"' for r in INSTALL_REQUIRES]
)

with open("$DOCKERFILE", "r") as fid:
  lines = fid.read().splitlines()

lines[-1] = f"RUN pip install --upgrade --upgrade-strategy only-if-needed {requires}"

with open("$DOCKERFILE", "w") as fid:
  fid.write( os.linesep.join(lines) )

EOF

docker build -t wbgt .

docker run -it \
      --rm \
      -p 8080:8888 \
      -p 8081:8787 \
      -v "${PWD}"/:/home/jovyan/work \
      -v $HOME/.aws/credentials:/home/jovyan/.aws/credentials:ro \
      wbgt 
