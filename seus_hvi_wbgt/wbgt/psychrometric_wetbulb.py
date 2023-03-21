from numpy import arctan
from .calc import relative_humidity

def stull( Ta, Td ):
  """
  Wet bulb temperature from Stull method

  This formula for wet bulb temperature appears in:
    Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and 
      Air Temperature, Journal of Applied Meteorology and Climatology, 
      50(11), 2267-2269. Retrieved Jul 20, 2022, from 
      https://journals.ametsoc.org/view/journals/apme/50/11/jamc-d-11-0143.1.xml

  Inputs:
    Ta (ndarray) : Ambient (dry bulb) temperature (degree C)
    Td (ndarray) : Dew point temperature (degree C)

  """

  RH = relative_humidity( Ta, Td ) * 100.0
  return Ta*arctan( 0.151977*(RH + 8.313659)**(1.0/2.0) ) +\
         arctan( Ta + RH ) - arctan( RH - 1.676331 ) +\
         0.00391838*RH**(3.0/2.0)*arctan( 0.023101*RH ) -\
         4.686035
