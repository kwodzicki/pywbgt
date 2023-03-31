FROM jupyter/datascience-notebook:latest

CMD pip install -e ./work; start-notebook.sh

USER root

RUN apt-get update
RUN apt-get install -y libgeos-dev libomp-dev

USER jovyan

RUN pip install --upgrade --upgrade-strategy only-if-needed "numpy" "scipy" "pandas" "metpy" "pvlib"
