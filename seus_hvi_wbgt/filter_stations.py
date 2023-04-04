"""
Filter stations to only long records

Based on proposal document, stations must meet certain criteria to be
included in the study; minimum number of days per year reporing for
minimum number of consecutive years.

Functions in this module enable the filtering of stations based on
these criteria 

"""

from ncei.ISD import utils, reader

def check_consecutive( data, syear, eyear, min_days ):
    """

    Arguments:
        data (DataFrame) : Station data
        syear (int) : Starting year for consecutive check
        eyear (int) : Ending year for consecutive check
        min_days (int)  : Number of days in a give year that must have data
            for year to be considered

    Returns:
        tuple : First element is list of years for and second is
            number of consecutive years that meet criteria. Not that
            this list 'resets' if a year does not meet the number of
            days criteria.

    """

    consecutive = [0]
    years       = []
    for year in range(syear, eyear+1):
        years.append( year )
        # hourly data so divide by 24 for number of days of data
        ndays = sum(data['datetime'].dt.year == year) / 24
        if ndays > min_days:
            consecutive.append( consecutive[-1] + 1 )
        else:
            consecutive.append( 0 )

    return years, consecutive[1:]

def check_station( data, syear, eyear = None, ndays = 300, nyears = 10 ):
    """

    Arguments:
        data (DataFrame) : Data from station to check if meets time criteria
        syear (int) : Starting year to check
 
    Keyword arguments:
        ndays (int)  : Number of days in a give year that must have data
            for year to be considered
        nyears (int) : Threshold for how many consecutive years must meet
            the days per year threshold

    """

    if eyear is None:
        eyear = syear

    _, consecutive = check_consecutive( data, syear, eyear, ndays )

    return any( c >= nyears for c in consecutive )

def main(wban, syear, eyear):
    """
    Check given wban stations

    """

    variables = ['date', 'time', 'air temperature', 'dew point temperature']

    url_gen = utils.genRemoteList( range(syear, eyear+1), wban=wban )
    return reader.read( url_gen, variables = variables )
