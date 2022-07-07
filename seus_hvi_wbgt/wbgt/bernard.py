import numpy

from metpy.units import units

def atmosVaporPres( temp, a = 0.6107, b = 17.27, c = 237.3 ):
  """
  Compute saturation vapor pressure from temperature

  From Bernard (1999)

  Arguments:
  temp (ndarry) : Temperature at which to compute saturation vapor pressure.
    Units of degree Celsius

  """

  return a * numpy.exp( b*temp / (temp + c) ) 

def psychrometricWetBulb(Ta, ea=None, Td=None, RH=None ):
  """
  Calculate psychrometricWetBulb from other variables

  Arguments:
    Ta (ndarray) : ambient temperature in degree Celsius

  Keyword arguments:
    ea (ndarray) : ambient vapor pressure in hectoPascals
    Td (ndarray) : Dew point temperature in degree Celsius
    RH (ndarray) : Relative humidity in percent

  """

  if ea is None:
    if Td is not None:
      ea = atmosVaporPres( Td )
    elif RH is not None:
      ea = RH/100.0 * atmosVaporPres( Ta )
    else:
      raise Exception( 'Must input one of "ea", "RH", or "Td"' )
 
  return 0.376 + 5.79*ea+ (0.388 - 0.0465*ea) * Ta
  
def factorC( wnd ):
  """Compute factor for natural wet bulb temperature

  Arguments:
    wnd (ndarray) : Wind speed in meters/second

  """

  C   = numpy.where( wnd < 0.03, 0.85, 1.0 )                                    # Where wind is less than 0.03, set values to 0.85, else set values to 1.0
  return numpy.where( wnd > 3.0, C, 0.96 + 0.069*numpy.log10( wnd ) )           # Where wind > 3.0, keep values of C, else compute C and return values 

def factore( wnd ):

  e = numpy.where( wnd < 0.1, 1.1, -0.1 )                                       # Where wind is less than 0.1, set values to 1.1, else set values to -0.1 
  return numpy.where( wnd > 1.0, e, 0.1/wnd**1.1 - 0.2 )                        # Where wind > 1.0, keep values of e, else compute e and return values 

def naturalWetBulb( Ta, Tpwb, Tg, wnd ):
  """
  Compute natural wet bulb temperature using Bernard method

  Method obtained from:
    https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GH000231

  """

  return numpy.where( (Tg-Ta) < 4.0,
    Ta - factorC( wnd )*(Ta - Tpwb),
    Tpwb + 0.25*(Tg-Ta) + factore( wnd )
  )

def bernard( u, Ta, Td ):
  """
  Compute WBGT using Bernard Method

  Arguments:
    u (Quantity) : Wind speed
    Ta (Qunatity) : Ambient temperature
    Td (Quantity) : Dew point temperature

  Returns:
    Quantity : Wet bulb globe temperature

  """

  u    = u.to('m/s')
  Ta   = Ta.to('degree_Celsius')
  Td   = Td.to('degree_Celsius')

  Tpwb = psychrometricWetBulb( Ta, Td = Td )
  Td   = globeTemperature( )

  return units.Quantity(
    0.7 * naturalWetBulb( Ta, Tpwb, Tg, u ) +\
    0.2 * Tg  +\
    0.1 * Ta,
    'degree_Celsius'
  )
