"""
Wet-bulb Globe Temperature using the Ono and Tonouchi (2012) method

"""

import numpy as np

from .calc import relative_humidity

def ono( datetime, lat, lon,
        solar, pres, temp_air, temp_dew, speed,
        f_db=None,
        cosz=None,
        **kwargs,
    ):
    """
    Compute WBGT using Dimiceli method

    Arguments:
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        lat (float) : Latitude of observations
        lon (float) : Longitude of observations
        solar (Quantity) : Solar irradiance; unit of power over area
        pres (Quantity) : Atmospheric pressure; unit of pressure
        temp_air (Quantity) : Ambient temperature; unit of temperature
        temp_dew (Quantity) : Dew point temperature; unit of temperature
        speed (Quantity) : Wind speed; units of speed

    Keyword arguments:
        f_db (float) : Direct beam radiation from the sun; fraction
        cosz (float) : Cosine of solar zenith angle
        wetbulb (str) : Name of wet bulb algorithm to use:
            {dimiceli, stull} DEFAULT = dimiceli
        natural_wetbulb (str) : Name of the natural wet bulb algorithm to use:
            {malchaire, hunter_minyard} DEFAULT = malchaire

    Returns:
        dict :
            - Tg : Globe temperatures in ndarray
            - Tpsy : psychrometric wet bulb temperatures in ndarray
            - Tnwb : Natural wet bulb temperatures in ndarray
            - Twbg : Wet bulb-globe temperatures in ndarray

    """

    solar    =    solar.to( 'kilowatt/m**2'     ).magnitude
    pres     =     pres.to( 'hPa'               ).magnitude
    temp_air = temp_air.to( 'degree_Celsius'    ).magnitude
    temp_dew = temp_dew.to( 'degree_Celsius'    ).magnitude
    speed    =    speed.to( 'meters per second' ).magnitude

    relhum = relative_humidity( temp_air, temp_dew ) * 100.0

    wbgt = (
        0.73500 * temp_air +
        0.03740 * relhum +
        0.00292 * temp_air*relhum +
        7.61900 * solar -
        4.55700 * solar**2 -
        0.05720 * speed -
        4.064
    )

    return {
        'Tg'   : np.nan,
        'Tpsy' : np.nan,
        'Tnwb' : np.nan,
        'Twbg' : wbgt
    }
