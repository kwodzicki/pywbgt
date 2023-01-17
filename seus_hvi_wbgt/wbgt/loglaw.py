from numpy import log

def loglaw(vRef, zRef, zNew=2.0, z0=0.4):
  """
  To downscale wind using Log Law

  Arguments:
    vRef (pint.Quanity) : Known velocity at height zRef with units
      of velocity
    zRef (pint.Quanityt) : Reference height where vRef is known with 
      units of distance

  Keyword arguments:
    zNew (float) : New height (in meters) above ground level
    Z0 (ndarray) : Roughness length in the current wind direction

  Returns:
    ndarray : Velocity calculated at height zNew

  """

  return (
    vRef.to('meters per second') 
    * log( zNew  / z0 ) 
    / log( zRef.to('meters').magnitude  / z0 )
  )

