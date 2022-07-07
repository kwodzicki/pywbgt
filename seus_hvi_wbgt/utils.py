from numpy import deg2rad, exp, cos

from pyorbital.astronomy import sun_zenith_angle, sun_earth_distance_correction

def f_dir_z( date, lat, lon, S, S0 = 1367.0 ):
  """
  Compute direct beam solar fraction and zenith angle

  Arguments:
    date (datetime) : A timezone aware datetime object
    lat (float) : Latitude of location; positive to north
    lon (float) : Longitude of location; positive to east
    S (float) : Global solar radiation in W/m**2

  Keyword arguments:
    S0 (float) : Solar constant; should not need changed

  Returns:
    tuple : 2 floats that are the Fraction of global solar 
      radiation that is direct beam and solar zenith angle,
      respectively

  """

  theta = sun_zenith_angle( date, lon, lat )
  if theta > 89.5:
    return 0.0

  theta = deg2rad( theta )
  d     = sun_earth_distance_correction( date )
  Ss    = S / (S0 * cos( theta ) / d**2)

  return exp( 3 - 1.34*Ss - 1.65/Ss ), theta
  
