"""
Utilities used across modules

Common utilities used across the pywbgt packge are found here

"""

from typing import Union

import numpy as np
from pandas import to_datetime, to_timedelta, DatetimeIndex
from xarray import DataArray


def datetime_adjust(datetime, gmt, avg):
    """
    Adjust datetime based on averaging window and gmt

    Arguments:
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        gmt (ndarray) LST-GMT difference  (hours; negative in USA)
        avg (ndarray) : averaging time of the meteorological inputs (minutes)

    Returns:
        pandas.DatetimeIndex : Adjusted datetimes based on the avg and
            gmt values

    """

    datetime = datetime_check(datetime)
    # Set default data averaging interval and compute time offset so that
    # time is in the middle of the sampling interval
    if avg is None:
        avg = 1.0
    dt = to_timedelta(avg / 2.0, 'minute')

    # If gmt is NOT None (i.e., it is set), then adjust the time delta
    if gmt is not None:
        dt = dt + to_timedelta(gmt, 'hour')

    # Adjust time using time delta
    return datetime - dt


def datetime_check(datetime):
    """
    Attempt to convert datetime to a DatetimeIndex object

    """

    if isinstance(datetime, DataArray):
        return datetime.values

    if isinstance(datetime, DatetimeIndex):
        return datetime

    if np.issubdtype(datetime.dtype, np.datetime64):
        return datetime

    return to_datetime(datetime)


def convert_units(val: Union['DataArray', 'ndarray'], unit_str: str):

    # If is an xarray object with metpy attribute
    if hasattr(val, 'metpy'):
        return (
            val
            .metpy
            .convert_units(unit_str)
            .metpy
            .dequantify()
        )

    return val.to(unit_str).magnitude
