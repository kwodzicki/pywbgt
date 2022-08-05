from numpy import exp

def saturation_vapor_pressure( T ):
  """
  Compute saturation vapor pressure

  Uses the Bolton (1980) formula to calculate saturation vapor
  pressure given a temperature in degrees Celsius

  Arguments:
    T (ndarray) : Temperature in degrees Celsius

  Returns:
    ndarray : Staturation vapor pressure(s) in hPa

  """

  return 6.112 * exp( 17.67 * T / (T + 243.5) )
 
def relative_humidity( Ta, Td ):
  """
  Compute relative humidity given temperature and dew point

  The relative humidity is computed by dividing the 
  saturation vapor pressure at dew point temperature by
  the saturation vapor pressure at the air temperature.

  Arguments:
    Ta (ndarray) : Ambient temperature in degrees Celsius
    Td (ndarray) : Dew point temperature in degrees Celsius

  Returns:
    ndarray : Relative humidity as dimensionless fraction

  """
 
  return saturation_vapor_pressure( Td )/saturation_vapor_pressure( Ta )
