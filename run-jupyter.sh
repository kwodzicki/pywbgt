#!/bin/bash

docker build -t wbgt .

#      -p 8081:8787 \

docker run -it \
      --rm \
      -p 8080:8888 \
      -p 8081:8787 \
      -v "${PWD}"/:/home/jovyan/work \
      -v $HOME/.aws/credentials:/home/jovyan/.aws/credentials:ro \
      wbgt 
