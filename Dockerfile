FROM jupyter/minimal-notebook:python-3.10

CMD pip install -e ./work; start-notebook.sh

USER root

RUN apt-get update
RUN apt-get install -y gcc g++ libgeos-dev libomp-dev
RUN apt-get install -y dvipng texlive-latex-extra texlive-fonts-recommended cm-super

USER jovyan

RUN pip install --upgrade --upgrade-strategy only-if-needed "cython" "numpy" "scipy" "pandas" "pyarrow" "metpy" "pvlib" "matplotlib" 