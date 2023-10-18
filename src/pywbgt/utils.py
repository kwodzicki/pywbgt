from pandas import to_datetime, DatetimeIndex

def datetime_check(datetime):
    """
    Attempt to convert datetime to a DatetimeIndex object

    """

    if isinstance(datetime, DatetimeIndex):
        return datetime

    return to_datetime(datetime)
 
