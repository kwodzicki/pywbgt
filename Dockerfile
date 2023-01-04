FROM jupyter/datascience-notebook:latest

# Install Rust
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

RUN conda install cartopy
RUN pip install zarr pyarrow fastparquet \
xarray pandas \
matplotlib polars \
boto3 s3fs \
pvlib \
metpy

USER root
RUN apt install libgomp1

USER jovyan
ADD --chown=jovyan:users . src
RUN pip install src/
