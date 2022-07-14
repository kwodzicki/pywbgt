import numpy

from .liljegren import solar_parameters

#from ..utils import f_dir_z
from .. import SIGMA

#from .utils import wetBulb

def vaporPressure( T ):

  return 6.112 * numpy.exp( 17.67 * T / (T + 243.5) )
 
def relativeHumidity( Ta, Td ):

  return vaporPressure( Td ) / vaporPressure( Ta ) * 100.0

def wetBulbRH( Ta, RH ): 

  return -5.806    + 0.672   *Ta -  0.006   *Ta**2       +\
       (  0.061    + 0.004   *Ta + 99.000e-6*Ta**2) * RH +\
       (-33.000e-6 - 5.000e-6*Ta -  1.000e-7*Ta**2) * RH**2
 
def wetBulb( Ta, Td ): 

  return wetBulbRH( Ta, relativeHumidity( Ta, Td ) )
   
def chfc( S=None, Z=None):
  """
  Convective heat flow coefficient

  Returns:
    float : THe convective heat flow coefficient

  Notes:
    Default value obtained from:
    https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf

  """

  if (S is not None) and (Z is not None):
    return a * S**b * numpy.cos( Z )**c

  return 0.315

def atmosVaporPres( Ta, Td, P ):
  """
  Compute atmospheric vapor pressure

    Ta (float) : ambient temperature in degrees Celsius
    Td (float) : dew point temperature in degrees Celsius
    P (float) : Barometric pressure in hPa

  Returns:
    float : atmospheric vapor pressure

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

  return 0.575 * atmosVaporPres( Ta, Td, P )**(1.0/7.0)


def factorB( Ta, Td, P, S, f_db, f_dif, z ):
  """

  Arguments:

  """

  return S * ( f_db / ( 4.0*SIGMA*numpy.cos(z) ) + 1.2 * f_dif / SIGMA ) +\
      thermalEmissivity( Ta, Td, P ) * Ta**4

def factorC( u ):
  """

  Arguments:
    u (float) : wind speed in meters per hour

  """

  return chfc() * u**0.58 / 5.3865e-8

def globeTemperature( Ta, Td, u, P, S, f_db, f_dif, z ):
  """
  Compute globe temperature

  Arguments:
    Ta (float) : ambient temperature in degrees Celsius
    Td (float) : dew point temperature in degrees Celsius
    u (float) : wind speed in meters per hour
    P (float) : Barometric pressure in inHg 
    S (float) : solar irradiance in Watts per meter**2
    f_dif (float) : diffuse radiation from the sun
    Z (float) : solar zenith angle

  Notes:
    Chapter 26 of IAENG Transactions on Engineering Technologies: 
      "Black Globe Temperature Estimate for the WBGT Index"

    https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf

  """

  B = factorB( Ta, Td, P, S, f_db, f_dif, z)
  C = factorC( u )

  return (B + C*Ta + 7.68e6) / (C + 2.56e5)
 
def dimiceli( lat, lon, 
              year, month, day, hour, minute, 
              solar, pres, Tair, Tdew, speed, f_db=None, z=None, **kwargs ):
  """
  Compute WBGT using Dimiceli method

  Arguments:
    lat (float) : Latitude of observations
    lon (float) : Longitude of observations
    year (ndarray) : UTC Year of the data values
    month (ndarray) : UTC Month of the data values
    day (ndarray) : UTC Day of the data values
    hour (ndarray) : UTC Hour of the data values; can be any time zone as long
    minute (ndarray) : UTC Minute of the data values
    u (Quantity) : Wind speed; units of speed
    Ta (Quantity) : Ambient temperature; unit of temperature
    Td (Quantity) : Dew point temperature; unit of temperature
    P (Quantity) : Atmospheric pressure; unit of pressure
    S (Quantity) : Solar irradiance; unit of power over area

  Keyword arguments:
    f_db (float) : Direct beam radiation from the sun; fraction
    z (float) : Solar zenith angle in radians

  Returns:
    tuple : Globe, natural wet bulb, and wet bulb globe temperatures

  """

  solar = solar.to( 'watt/m**2'       ).magnitude
  pres  =  pres.to( 'inHg'            ).magnitude
  Tair  =  Tair.to( 'degree_Celsius'  ).magnitude
  Tdew  =  Tdew.to( 'degree_Celsius'  ).magnitude
  speed = speed.to( 'meters per hour' ).magnitude

  if (f_db is None) or (z is None):
    solar = solar_parameters( lat, lon, 
                              year, month, day, hour, minute, 
                              solar )
    if f_db is None: f_db = solar[2]
    if z    is None: z    = numpy.arccos( solar[1] ) 
    solar = solar[0]

  f_dif = 1.0 - f_db
  Tg    = globeTemperature( Tair, Tdew, speed, pres, solar, f_db, f_dif, z )
  Tnwb  = wetBulb( Tair, Tdew )

  return Tg, Tnwb, 0.7*Tnwb + 0.2*Tg + 0.1*Tair
