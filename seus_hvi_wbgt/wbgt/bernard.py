import numpy
from scipy.optimize import fsolve
from metpy.calc import relative_humidity_from_dewpoint as rhTd

from . import liljegren 
from . import SIGMA, EPSILON

"""
References

Bernard, T. E., & Cross, R. R. (1999). 
    Heat stress management: Case study in an aluminum smelter. 
    International journal of industrial ergonomics, 23(5-6), 609-620.

Bernard, T.E., & Pourmoghani, M. (1999). 
    Prediction of workplace wet bulb global temperature. 
    Applied occupational and environmental hygiene, 14 2, 126-34 .

"""

def chtc( Tg, Ta, u ):
  """
  Convective heat transfer coefficent

  This is equation 5 from Bernard and Pourmoghani (1999) used to
  compute the convective heat transfer coefficient.

  Arguments:
    Tg (ndarray) : Globe temperature; degree Celsius
    Ta (ndarray) : Ambient temperature; degree Celsius
    u (ndarray) : Wind speed; meters/second

  Kewyord arguments:
    None.

  Returns:
    ndarray : Convective heat transfer coefficient values

  """

  return ( (10.9*u**0.566)**3 + (0.35 + 1.77*(Tg-Ta)**0.25)**3 )**(1.0/3.0)

def funcTg( Tg, Ta, u, P, S, epsilon ):#f_db, f_dif, z )
  """
  Differential equation for globe temperature

  This is equation 4 from Bernard and Pourmoghani (1999) used to 
  determine the globe temperature through an iterative solver.

  Arguments:
    Tg (ndarray) : First guess for globe temperature; degree Celsius
    Ta (ndarray) : Ambient temperature; degree Celsius
    u (ndarray) : Wind speed; meters/second
    P (ndarray) : Atmospheric pressure; hPa
    S (ndarray) : Radiant heat flux incident on the globe;
      currently assuming this to be the solar irradiance; W/m**2
    epsilon (ndarray) : Emissivity of the the globe; default value
      will be same as that used for the Liljegren method.

  """

  hg = chtc( Tg, Ta, u )                                                        # Compute the convective heat transfer coefficient for given values
  return hg*(Tg-Ta) - epsilon*(S - SIGMA*Tg**4)                                 # Return result of the differential equation for determining globe temperature

def factorC( wnd ):
  """
  Compute factor for natural wet bulb temperature without radiant heat

  This factor is used to related psychrometric wet bulb temperature
  to natural wet bulb temperature in the absence of radiant heat.

  Arguments:
    wnd (ndarray) : Wind speed in meters/second

  """

  C   = numpy.where( wnd < 0.03, 0.85, 0.96 + 0.069*numpy.log10( wnd ) )                                    # Where wind is less than 0.03, set values to 0.85, else set values to 1.0
  return numpy.where( wnd > 3.0, 1.0, C )           # Where wind > 3.0, keep values of C, else compute C and return values 


def factore( wnd ):
  """
  Compute factor for natural wet bulb temperature with radiant heat

  This factor is used to related psychrometric wet bulb temperature
  to natural wet bulb temperature in the presence o f radiant heat.

  Arguments:
    wnd (ndarray) : Wind speed in meters/second

  """

  e = numpy.where( wnd < 0.1, 1.1, 0.1/wnd**1.1 - 0.2 )                                       # Where wind is less than 0.1, set values to 1.1, else set values to -0.1 
  return numpy.where( wnd > 1.0, -0.1, e )                        # Where wind > 1.0, keep values of e, else compute e and return values 

def saturation_vapor_pressure( T ):
  """
  Compute saturation vapor pressure from temperature

  Equation taken from From Bernard and Pourmoghani (1999), with reference
  to Beranrd and Cross (1999)

  Arguments:
    T (ndarry) : Temperature at which to compute saturation vapor pressure.
      Units of degree Celsius

  Returns:
    ndarray : vapor pressure in kPa

  """

  return 0.6107*numpy.exp( 17.27*T / (T + 237.3) ) 

def globeTemperature( Ta, u, P, S, epsilon=EPSILON  ):
  """
  Determine globe temperature through iterative solver

  Aruguments:
    Ta (ndarray) : Ambient temperature; degree Celsius
    u (ndarray) : Wind speed; meters/second
    P (ndarray) : Atmospheric pressure; hPa
    S (ndarray) : Radiant heat flux incident on the globe;
      currently assuming this to be the solar irradiance; W/m**2

  Keyword arguments:
    epsilon (ndarray) : Emissivity of the the globe; default value
      will be same as that used for the Liljegren method.

  Returns:
    ndarray : Globe temperature determined through iterative solver

  """

  return fsolve( funcTg, Ta+1.0, args=(Ta, u, P, S, epsilon) )

def psychrometricWetBulb(Ta, ea=None, Td=None, RH=None ):
  """
  Calculate psychrometric wet bulb from other variables

  Arguments:
    Ta (ndarray) : ambient temperature in degree Celsius

  Keyword arguments:
    ea (ndarray) : ambient vapor pressure in kiloPascals
    Td (ndarray) : Dew point temperature in degree Celsius
    RH (ndarray) : Relative humidity in percent

  """

  if ea is None:
    if Td is not None:
      ea = saturation_vapor_pressure( Td )
    elif RH is not None:
      ea = RH/100.0 * saturation_vapor_pressure( Ta )
    else:
      raise Exception( 'Must input one of "ea", "RH", or "Td"' )
 
  return 0.376 + 5.79*ea + (0.388 - 0.0465*ea)*Ta
  

def naturalWetBulb( Ta, Tpsy, Tg, wnd ):
  """
  Compute natural wet bulb temperature using Bernard method

  If the globe temperature exceeds the dry bulb temperature by more than 
  4 degree C, then the effect of radiant heat is considered in the 
  calculation (factor e is computed). If this difference is 4 degree C or
  less, then no radiant heat effect is considered (factor C is computed)

  Arguments:
    Ta (ndarray) : Ambient temperature; degree Celsius
    Tpsy (ndarray) : Psychrometric wet bulb temperature; degree Celsius
    Tg (ndarray) : Globe temperature; degree Celsius
    wnd (ndarray) : Wind speed; meters/second

  Keyword arguments:
    None.

  Returns:
    ndarray : Natural wet bulb temperature in degree Celsius

  Method obtained from:
    https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GH000231

  Notes:
    In the reference paper, the psychrometric wet bulb temperture is 
    Tpwb, however, to be consistent with other papers, and because I
    believe it makes more sense, we use Tpsy for psychrometric wet bulb

  """

  return numpy.where( (Tg-Ta) < 4.0,
    Ta - factorC( wnd )*(Ta - Tpsy),
    Tpsy + 0.25*(Tg-Ta) + factore( wnd )
  )

def bernard( lat, lon,
             year, month, day, hour, minute,
             solar, pres, Tair, Tdew, speed, **kwargs):
  """
  Compute WBGT using Bernard Method

  Arguments:
    lat (float) : Latitude of observations
    lon (float) : Longitude of observations
    year (ndarray) : UTC Year of the data values
    month (ndarray) : UTC Month of the data values
    day (ndarray) : UTC Day of the data values
    hour (ndarray) : UTC Hour of the data values; can be any time zone as long
    minute (ndarray) : UTC Minute of the data values
    solar (Quantity) : Solar irradiance; unit of power over area
    pres (Quantity) : Atmospheric pressure; unit of pressure
    Tair (Quantity) : Ambient temperature; unit of temperature
    Tdew (Quantity) : Dew point temperature; unit of temperature
    speed (Quantity) : Wind speed; units of speed

  Returns:
    tuple: Globe, natural wet bulb, psychrometric wet bulb, and
      wet bulb globe temperatures

  """

  solar, cza, fdir = liljegren.solar_parameters(
    lat, lon, year, month, day, hour, minute, solar.to('watt/m**2').magnitude
  )

  pres   = pres.to('hPa').magnitude
  speed  = speed.to('meter/second').magnitude
  relhum = rhTd( Tair, Tdew ).to('percent').magnitude

  Tg   = liljegren.globeTemperature( 
           Tair.to('kelvin').magnitude, relhum/100.0, 
           pres, speed, solar, fdir, cza 
  )
  Tg[ Tg <= -9990.0] = float('nan')
  Tair = Tair.to('degree_Celsius').magnitude
  Tpsy = psychrometricWetBulb( Tair, RH = relhum )
  Tnwb = naturalWetBulb( Tair, Tpsy, Tg, speed )

  return {
   'Tg'   : Tg, 
   'Tnwb' : Tnwb, 
   'Tpsy' : Tpsy, 
   'Twbg' : 0.7*Tnwb + 0.2*Tg + 0.1*Tair
  }
