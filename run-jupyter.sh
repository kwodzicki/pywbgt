#!/bin/bash

python3 ./update_dockerfile.py ./Dockerfile

docker build -t wbgt .

docker run -it \
    --rm \
    -p 8080:8888 \
    -p 8081:8787 \
    -v "${PWD}"/:/home/jovyan/work \
    -v $HOME/.aws/credentials:/home/jovyan/.aws/credentials:ro \
    "$@" \
    wbgt 
