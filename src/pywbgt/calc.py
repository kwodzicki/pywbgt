"""
Calculate various atmospheric properties

"""

from numpy import exp, log

def loglaw(velo_ref, z_ref, z_new=2.0, z_rough=0.1, zp_displace=0.0):
    """
    To downscale wind using Log Law

    Arguments:
        velo_ref (pint.Quanity) : Known velocity at height z_ref with units
            of velocity
        z_ref (pint.Quanityt) : Reference height where velo_Ref is known with 
            units of distance

    Keyword arguments:
        z_new (float) : New height (in meters) above ground level
        z_rough (float,ndarray) : Roughness length (in meters) in the current
            wind direction
        zp_displace (float, ndarray) : Zero-plane displacement;
            height in meters

    Returns:
        ndarray : Velocity calculated at height z_new

    """

    return (
      velo_ref.to('meters per second')
      * log( (z_new-zp_displace)  / z_rough )
      / log( (z_ref.to('meters').magnitude-zp_displace)  / z_rough )
    )
