"""
Calculate various atmospheric properties

"""

from numpy import exp, log

def saturation_vapor_pressure( temperature ):
    """
    Compute saturation vapor pressure
  
    Uses the Bolton (1980) formula to calculate saturation vapor
    pressure given a temperature in degrees Celsius
  
    Arguments:
        temperature (ndarray) : Temperature in degrees Celsius
  
    Returns:
        ndarray : Staturation vapor pressure(s) in hPa
  
    """

    return 6.112 * exp( 17.67 * temperature / (temperature + 243.5) )

def relative_humidity( temp_a, temp_d ):
    """
    Compute relative humidity given temperature and dew point
  
    The relative humidity is computed by dividing the 
    saturation vapor pressure at dew point temperature by
    the saturation vapor pressure at the air temperature.

    Arguments:
        temp_a (ndarray) : Ambient temperature in degrees Celsius
        temp_d (ndarray) : Dew point temperature in degrees Celsius

    Returns:
        ndarray : Relative humidity as dimensionless fraction

    """

    return saturation_vapor_pressure( temp_d )/saturation_vapor_pressure( temp_a )

def loglaw(velo_ref, z_ref, z_new=2.0, z_rough=0.4):
    """
    To downscale wind using Log Law

    Arguments:
        velo_ref (pint.Quanity) : Known velocity at height z_ref with units
            of velocity
        z_ref (pint.Quanityt) : Reference height where velo_Ref is known with 
            units of distance

    Keyword arguments:
        z_new (float) : New height (in meters) above ground level
        z_rough (ndarray) : Roughness length in the current wind direction

    Returns:
        ndarray : Velocity calculated at height z_new

    """

    return (
      velo_ref.to('meters per second')
      * log( z_new  / z_rough )
      / log( z_ref.to('meters').magnitude  / z_rough )
    )
