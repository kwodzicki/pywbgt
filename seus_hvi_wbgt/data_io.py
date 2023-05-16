"""
General IO utilities

Adds reader/writer for parquet files that writes
variable metadata to the parquet files.

"""

import logging

from pandas import read_parquet as pd_read_parquet
from pyarrow import parquet, Schema

def read_parquet( file, **kwargs ):
    """
    Read data from parquet file

    Data from a parquet file is read into a pandas
    DataFrame, including attributes for variables
    (columns).

    Arguments:
        file (str) : Full-path of parquet file to read in

    Keywords arguments:
        **kwargs : passed to pandas.read_parquet

    Returns:
        tuple : pandas.DataFrame containing data and
            dict containing variable (column) attributes.

    """

    data   = pd_read_parquet(file, **kwargs)
    schema = parquet.read_schema( file )
    attrs  = {}
    for fname in schema.names:
        metadata = schema.field(fname).metadata
        if not isinstance(metadata, dict):
            continue
        attrs[fname] = {
            key.decode() : val.decode()
            for key, val in metadata.items()
        }

    return data, attrs

def write_parquet( file, dataframe, attrs, **kwargs ):
    """
    Write DataFrame to parquet file

    A pandas DataFrame is written to a parquet file, including attributes
    for variables (columns).

    Arguments:
        file (str) : Full-path of parquet file to write data to
        dataframe (pandas.DataFrame) : Data to write to file
        attrs (dict) : Variable (column) attributes to write to the
            parquet file. Layout of dictionary is keys that match variable
            (column) names, with values being dictionaries of attribute 
            names (keys) and attribute values (value).
        **kwargs : Passed to pandas.DataFrame.to_parquet(). Note that
            the 'engine' keyword has not effect as it is forced
            to 'pyarrow'

    Returns:
        None.

    """

    log = logging.getLogger(__name__)
    _   = kwargs.pop('engine', None)

    schema = Schema.from_pandas( dataframe )
    for key, val in attrs.items():
        idx = schema.get_field_index( key )
        if idx == -1:
            log.error( "Failed to get schema field for '%s'", key )
            continue

        schema = schema.set(
            idx,
            schema.field(idx).with_metadata( val ),
        )

    dataframe.to_parquet(
        file,
        schema = schema,
        engine = 'pyarrow',
        **kwargs,
    )
