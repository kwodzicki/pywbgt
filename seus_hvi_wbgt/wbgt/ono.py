"""
Wet-bulb Globe Temperature using the Ono and Tonouchi (2012) method

"""

import numpy as np

from .utils import relative_humidity

def ono( lat, lon, datetime,
              solar, pres, Tair, Tdew, speed, f_db=None, cosz=None, **kwargs ):
  """
  Compute WBGT using Dimiceli method

  Arguments:
    lat (float) : Latitude of observations
    lon (float) : Longitude of observations
    datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
    solar (Quantity) : Solar irradiance; unit of power over area
    pres (Quantity) : Atmospheric pressure; unit of pressure
    Tair (Quantity) : Ambient temperature; unit of temperature
    Tdew (Quantity) : Dew point temperature; unit of temperature
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

  solar = solar.to( 'kilowatt/m**2'     ).magnitude
  pres  =  pres.to( 'hPa'               ).magnitude
  Tair  =  Tair.to( 'degree_Celsius'    ).magnitude
  Tdew  =  Tdew.to( 'degree_Celsius'    ).magnitude
  speed = speed.to( 'meters per second' ).magnitude

  rh = relative_humidity( Tair, Tdew ) * 100.0

  wbgt = (
    0.735*Tair +
    0.0374*rh + 
    0.00292*Tair*rh +
    7.619*solar - 
    4.557*solar**2 -
    0.0572*speed -
    4.064
  )

  return {
    'Tg'   : np.nan,
    'Tpsy' : np.nan,
    'Tnwb' : np.nan,
    'Twbg' : wbgt
  }

