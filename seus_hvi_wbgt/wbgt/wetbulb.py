from numpy import arctan

from .utils import relative_humidity

def dimiceli( Ta, RH ): 
  """
  Wet bulb temperature from Dimiceli method

  This formula for wet bulb temperature appears at the bottom of 
  "Estimation of Black Globe Temperature for Calculation of the WBGT Index"
  by Dimiceli and Piltz.

  https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf

  Inputs:
    Ta (ndarray) : Ambient (dry bulb) temperature (degree C)
    RH (ndarray) : Relative humidity (%)

  """

  return -5.806    + 0.672   *Ta -  0.006   *Ta**2       +\
       (  0.061    + 0.004   *Ta + 99.000e-6*Ta**2) * RH +\
       (-33.000e-6 - 5.000e-6*Ta -  1.000e-7*Ta**2) * RH**2

def stull( Ta, RH ):
  """
  Wet bulb temperature from Stull method

  This formula for wet bulb temperature appears in:
    Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and 
      Air Temperature, Journal of Applied Meteorology and Climatology, 
      50(11), 2267-2269. Retrieved Jul 20, 2022, from 
      https://journals.ametsoc.org/view/journals/apme/50/11/jamc-d-11-0143.1.xml

  Inputs:
    Ta (ndarray) : Ambient (dry bulb) temperature (degree C)
    RH (ndarray) : Relative humidity (%)

  """

  return Ta*arctan( 0.151977*(RH + 8.313659)**(1.0/2.0) ) +\
         arctan( Ta + RH ) - arctan( RH - 1.676331 ) +\
         0.00391838*RH**(3.0/2.0)*arctan( 0.023101*RH ) -\
         4.686035

def wetBulb( Ta, Td, method='stull' ): 
  """
  Compute wet bulb temperature from temperature and dew point 

  Two methods of computing the wet bulb temperature are possible using
  this function; a method from Stull (2011) and a method from Dimiceli.
  By defualt use the Stull method as that seems more robust/widely 
  accepted and is presented as the formula used on a poster by Sean Heuser[1]

  Inputs:
    Ta (ndarray) : Ambient (dry bulb) temperature (degree C)
    Td (ndarray) : Dew point temperature (degree C)

  Keyword arguments:
    method (str) : Set to either 'stull' or 'dimiceli' to choose method
      for wet bulb temperature calcualation. Case insensitive.

  References:
    [1] Heuser, S. (?). Verification of Wet Bulb Globe Temperature Index
      for the North Carolina ECONet. NCSU, Retreived Jul 20, 2022, from 
      https://legacy.climate.ncsu.edu/downloads/econet/WBGT.pdf

  """

  method = method.upper()
  if method == 'STULL':
    return stull( Ta, relative_humidity( Ta, Td )*100.0 )
  elif method == 'DIMICELI':
    return dimiceli( Ta, relative_humidity( Ta, Td )*100.0 )
  
  raise Exception( f'Unrecognized method : {method}!' )
