FROM jupyter/datascience-notebook:latest

# Install Rust
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install Polars
RUN conda install dask distributed \
zarr pyarrow fastparquet \
xarray pandas \
matplotlib polars cartopy \
awswrangler boto3 s3fs \
metpy

run pip install ./setup.py
