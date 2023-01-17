from .utils import relative_humidity

def malchaire(Ta, Td, Tw, Tg ):
  """
  Compute natural wet bulb temperature

  An algorithm developed by Malchaire (1976) is used to compute the
  natural wet bulb temperature using the wet bulb, globe, ambient, and
  dew point temperature. Relative humidty is computerd using the
  ambient and dew point temperatures.

  Arguments:
    Ta (ndarray) : Ambient temp (C)
    Td (ndarray) : Dew point temperature (C)
    Tw (ndarray) : Wet bulb (C)
    Tg (ndarray) : Globe temp (C)

  Reference:
    Malchaire, J. B., EVALUATION OF NATURAL WET BULB AND WET GLOBE THERMOMETERS.
        The Annals of Occupational Hygiene, Volume 19, Issue 3-4, December 1976,
        Pages 251–258, https://doi.org/10.1093/annhyg/19.3-4.251

  """

  rh = relative_humidity( Ta, Td )
  return (0.16*(Tg-Ta) + 0.8)/200.0 * (560.0 - 2.0*rh - 5.0*Ta) - 0.8 + Tw

def hunter_minyard(Tw, S, u):
  """
  Compute natural wet bulb temperature

  An algorithm shown in Hunter and Minyard (1999) to compute the 
  natural wetbulb temperture using the psychrometric wetbulb, 
  solar irradiance, and wind speed.
 
  Arguments:
    Tw (ndarray) : Wet bulb globe temperature; degree C
    S (ndarray) : Solar irradiance; W/m**2
    u (ndarray) : Wind speed; m/s

  Reference:
    Hunter, Charles H., and C. Olivia Minyard. 
        Estimating wet bulb globe temperature using standard meteorological
        measurements." Proceedings of the Conference: 2nd Conference on
        Environmental Applications, Long Beach, CA, USA. Vol. 18. 1999.
 
  """

  #return Tw + 0.021*S - 0.42*u + 1.93
  return Tw + 0.0021*S - 0.43*u + 1.93 # Adjustment made based on the formula in the Boyer paper; see nws_boyer()

def nws_boyer( Ta, Tw, S, u ):
  """
  Compute natural wet bulb temperature

  An algorithm shown in National Weather Service (NWS) documention
  from Broyer to compute the natural wetbulb temperture using
  psychrometric wetbulb, webulb depression,
  solar irradiance, and wind speed.
 
  Arguments:
    Ta (ndarray) : Wet bulb globe temperature; degree C
    Tw (ndarray) : Wet bulb globe temperature; degree C
    S (ndarray) : Solar irradiance; W/m**2
    u (ndarray) : Wind speed; m/s

  Reference:
    Boyer, Timothy R.
      NDFD Wet Bulb Globe Temperature Algorithm and Software Design
      https://vlab.noaa.gov/documents/6609493/7858379/
      NDFD+WBGT+Description+Document.pdf/fb89cc3a-0536-111f-f124-e4c93c746ef7?t=1642792547129 

  """

  return Tw + 0.001651*S - 0.09555*u + 0.13235*(Ta-Tw) + 0.20249
