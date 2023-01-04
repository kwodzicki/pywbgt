"""
WBGT from the Dimiceli method

"""

import numpy

from . import SIGMA
from .utils import relative_humidity
from .liljegren import solar_parameters
from .stull import psychrometricWetBulb as stullWetBulb
from .natural_wetbulb import malchaire, hunter_minyard 

def chfc( S=None, Z=None):
  """
  Convective heat flow coefficient

  Returns:
    float : The convective heat flow coefficient

  Notes:
    Default value obtained from:
    https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf

  """

  if (S is not None) and (Z is not None):
    return a * S**b * numpy.cos( Z )**c

  return 0.315

def atmospheric_vapor_pressure( Ta, Td, P ):
  """
  Compute atmospheric vapor pressure

    Ta (float) : ambient temperature in degrees Celsius
    Td (float) : dew point temperature in degrees Celsius
    P (float) : Barometric pressure in hPa

  Returns:
    float : atmospheric vapor pressure

  Note:
    This is an odd function that I cannot find a lot of informaiton
    about. I was able to determe that the 2nd and 3rd lines of the function
    come from Buck, A. 1981: New Equations for computing vapor pressure
    and enhancement factor and is an improved equation for calculating
    vapor pressure in air, as opposed to pure water vapor.

    However, the first line of the equation does not make much sense.
    The base of the equation is from Bolton (1980), but the use of
    (Td-Ta) is odd. It is some kind of differential pressure factor.
    Some quick testing indicates that this formula will ALWAYS give
    a slighly lower pressure value than the Bolton (1980) formula; however,
    the values tend to be farily similar as long as there is not a large
    differenc between Ta and Td

  """

  return numpy.exp( (17.67 * (Td - Ta) ) / (Td + 243.5) ) *\
    (1.0007 + 3.46e-6 * P) *\
    6.112 * numpy.exp( 17.502 * Ta / (240.97 + Ta) )

def thermalEmissivity( Ta, Td, P ): 
  """
  Compute thermal emissivity from readily available NWS values

  Arguments:
    Ta (float) : ambient temperature in degrees Celsius
    Td (float) : dew point temperature in degrees Celsius
    P (float) : Barometric pressure in hPa

  Returns:
    float : thermal emissivity

  """

  return 0.575 * atmospheric_vapor_pressure( Ta, Td, P )**(1.0/7.0)

def factorB( Ta, Td, P, S, f_db, cosz, **kwargs ):
  """

  Arguments:
    Ta (float) : ambient temperature in degrees Celsius
    Td (float) : dew point temperature in degrees Celsius
    P (float) : Barometric pressure in hPa 
    S (float) : solar irradiance in Watts per meter**2
    f_db (float) : Fraction of direct beam radiation
    cosz (float) : Cosine of solar zenith angle

  """
  
  f_dif = 1.0 - f_db
  return S * ( f_db/(4.0*SIGMA*cosz) + 1.2*f_dif/SIGMA ) +\
      thermalEmissivity( Ta, Td, P ) * Ta**4

def factorC( u, CHFC = None, **kwargs ):
  """

  Arguments:
    u (float) : wind speed in meters per hour

  """

  if CHFC is None: CHFC = chfc()
  return CHFC * u**0.58 / 5.3865e-8

def globeTemperature( Ta, Td, P, u, S, f_db, cosz, **kwargs ):
  """
  Compute globe temperature

  Arguments:
    Ta (float) : ambient temperature in degrees Celsius
    Td (float) : dew point temperature in degrees Celsius
    P (float) : Barometric pressure in hPa 
    u (float) : wind speed in meters per hour
    S (float) : solar irradiance in Watts per meter**2
    f_db (float) : Fraction of direct beam radiation
    cosz (float) : Cosine of solar zenith angle

  Notes:
    Chapter 26 of IAENG Transactions on Engineering Technologies: 
      "Black Globe Temperature Estimate for the WBGT Index"

    https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf

  """

  B = factorB( Ta, Td, P, S, f_db, cosz)
  C = factorC( u, **kwargs )

  return (B + C*Ta + 7.68e6) / (C + 2.56e5)

def psychrometricWetBulb( Ta, Td ): 
  """
  Wet bulb temperature from Dimiceli method

  This formula for wet bulb temperature appears at the bottom of 
  "Estimation of Black Globe Temperature for Calculation of the WBGT Index"
  by Dimiceli and Piltz.

  https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf

  Inputs:
    Ta (ndarray) : Ambient (dry bulb) temperature (degree C)
    Td (ndarray) : Dew point temperature (degree C)

  """
  
  RH = relative_humidity( Ta, Td ) * 100.0
  return -5.806    + 0.672   *Ta -  0.006   *Ta**2       +\
       (  0.061    + 0.004   *Ta + 99.000e-6*Ta**2) * RH +\
       (-33.000e-6 - 5.000e-6*Ta -  1.000e-7*Ta**2) * RH**2 

def dimiceli( lat, lon, datetime, 
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

  solar  = solar.to( 'watt/m**2'       ).magnitude
  pres   =  pres.to( 'hPa'             ).magnitude
  Tair   =  Tair.to( 'degree_Celsius'  ).magnitude
  Tdew   =  Tdew.to( 'degree_Celsius'  ).magnitude
  speed1 = speed.to( 'meters per hour' ).magnitude
  speed1 = numpy.clip( speed1, 1690.0, None )                                   # Ensure wind speed is at least one (1) mile/hour

  if (f_db is None) or (cosz is None):
    solar = solar_parameters( lat, lon, datetime, solar, **kwargs )
    if cosz is None: cosz = solar[1] 
    if f_db is None: f_db = solar[2]
    solar = solar[0]

  Tg     = globeTemperature( Tair, Tdew, pres, speed1, solar, f_db, cosz, **kwargs)
  wbMeth = kwargs.get('wetbulb', 'DIMICELI').upper()
  if wbMeth == 'DIMICELI':
    Tpsy = psychrometricWetBulb( Tair, Tdew) 
  elif wbMeth == 'STULL':
    Tpsy = stullWetBulb( Tair, Tdew)
  else:
    raise Exception( f"Invalid option for 'wetbulb' : {wbMeth}" )

  nwb   = kwargs.get('natural_wetbulb', 'MALCHAIRE').upper()
  if nwb == 'MALCHAIRE':
    Tnwb  = malchaire( Tair, Tdew, Tpsy, Tg )
  elif nwb == 'HUNTER_MINYARD':
    Tnwb  = hunter_minyard( Tpsy, solar*f_db, speed.to('meter per second').magnitude )
  else:
    raise Exception( f"Unsupported value for 'natural_wetbulb' : {nwb}" )

  return {
    'Tg'   : Tg,
    'Tpsy' : Tpsy, 
    'Tnwb' : Tnwb, 
    'Twbg' : 0.7*Tnwb + 0.2*Tg + 0.1*Tair
  }
