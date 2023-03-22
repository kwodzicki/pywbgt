FROM jupyter/datascience-notebook:latest

USER root
RUN apt-get update
RUN apt-get install -y libgeos-dev libomp-dev
#RUN apt install libgomp1

USER jovyan
RUN pip install cartopy
RUN pip install --upgrade --upgrade-strategy only-if-needed "numpy" "pandas" "metpy" "pvlib"
CMD pip install -e ./work; start-notebook.sh