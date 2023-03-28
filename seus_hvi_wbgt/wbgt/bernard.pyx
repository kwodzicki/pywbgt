"""
WBGT using the Bernard method

References

Bernard, T. E., & Cross, R. R. (1999). 
    Heat stress management: Case study in an aluminum smelter. 
    International journal of industrial ergonomics, 23(5-6), 609-620.

Bernard, T.E., & Pourmoghani, M. (1999). 
    Prediction of workplace wet bulb global temperature. 
    Applied occupational and environmental hygiene, 14 2, 126-34 .

"""

from cython.parallel import prange
import numpy

from libc.math cimport fabs, pow
cimport numpy
cimport cython

from metpy.units import units

from . import SIGMA
from .liljegren import solar_parameters
from .calc import saturation_vapor_pressure, loglaw

cdef:
    int   MAX_ITER  = 50
    float CtoK      = 273.15
    float CONVERGE  = 1.0E-3
    float ALPHA_SFC = 0.00
    float NaN       = numpy.nan
    float SIGMAB    = SIGMA

@cython.cdivision(True)
cdef float _conv_heat_trans_coeff( float temp_g, float temp_a, float speed) nogil:
    """
    Convective heat transfer coefficient

    Compute the convective heat transfer coefficient
    from Bernard method

    Arguments:
        temp_g : Globe temperature; Kelvin
        temp_a : Air temperature; Kelvin
        speed : Wind speed; meter/second

    Returns:
        convective heat transfer coefficient

    """

    cdef float coeff, delta_t = temp_g-temp_a
    coeff = pow(
        pow(10.9*pow(speed, 0.566), 3) + 
        pow(0.35 + 1.77*pow(fabs(delta_t), 0.25), 3),
        1.0/3.0
    )
    if (delta_t < 0.0): 
        return -coeff
    return coeff

@cython.cdivision(True)
cdef float _globe_temperature( float temp_a, float speed, float pres, float solar, float f_db, float cosz ) nogil:
    """
   Determine globe temperature through iterative solver
  
    Aruguments:
        temp_a (ndarray) : Ambient temperature; Kelvin 
        speed (ndarray) : Wind speed; meters/second
        pres (ndarray) : Atmospheric pressure; hPa
        solar (ndarray) : Radiant heat flux incident on the globe;
            currently assuming this to be the solar irradiance; W/m**2
  
    Keyword arguments:
  
    Returns:
        ndarray : Globe temperature determined through iterative solver
  
    """

    cdef:
        int ii
        float h, temp_g_new, temp_g = temp_a

    for ii in range( MAX_ITER ):
        h = _conv_heat_trans_coeff( temp_g, temp_a, speed )
        temp_g_new = pow(
            temp_a**4 +
            solar/SIGMAB/2.0 * (1 + f_db*(1/2.0/cosz - 1.0) + ALPHA_SFC) -
            h/SIGMAB * (temp_g-temp_a),
            0.25
        )

        if fabs(temp_g_new-temp_g) < CONVERGE: 
            return temp_g_new
        temp_g = 0.9*temp_g + 0.1*temp_g_new

    return NaN

def conv_heat_trans_coeff( temp_g, temp_a, speed ):
    """
    Convective heat transfer coefficent
  
    This is equation 5 from Bernard and Pourmoghani (1999) used to
    compute the convective heat transfer coefficient.
  
    Arguments:
        temp_g (ndarray) : Globe temperature; degree Celsius
        temp_a (ndarray) : Ambient temperature; degree Celsius
        speed (ndarray) : Wind speed; meters/second
  
    Kewyord arguments:
        None.
  
    Returns:
        ndarray : Convective heat transfer coefficient values
  
    """

    delta_t = temp_g-temp_a
    coeff = (
        (10.9*speed**0.566)**3 + (0.35 + 1.77*abs(delta_t)**0.25)**3
    )**(1.0/3.0)

    return numpy.where(
        delta_t < 0,
        -coeff,
        coeff,
    )

def factor_c( speed ):
    """
    Compute factor for natural wet bulb temperature without radiant heat
  
    This factor is used to related psychrometric wet bulb temperature
    to natural wet bulb temperature in the absence of radiant heat.
  
    Arguments:
        speed (ndarray) : Wind speed in meters/second
  
    """

    fac_c      = numpy.full( speed.shape, 0.85 )
    idx        = numpy.where( speed>= 0.03 )
    fac_c[idx] = 0.96 + 0.069*numpy.log10( speed[idx] )
    # Where wind > 3.0, keep values of C, else compute C and return values
    return numpy.where( speed > 3.0, 1.0, fac_c )

def factor_e( speed ):
    """
    Compute factor for natural wet bulb temperature with radiant heat
  
    This factor is used to related psychrometric wet bulb temperature
    to natural wet bulb temperature in the presence o f radiant heat.
  
    Arguments:
        speed (ndarray) : Wind speed in meters/second
  
    """

    fac_e      = numpy.full( speed .shape, 1.1 )
    idx        = numpy.where( speed >= 0.1 )
    fac_e[idx] = 0.1/speed[idx]**1.1 - 0.2
    # Where wind > 1.0, keep values of e, else compute e and return values
    return numpy.where( speed > 1.0, -0.1, fac_e )

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def globe_temperature( temp_a, speed, pres, solar, f_db, cosz, ):
    """
    Determine globe temperature through iterative solver
  
    Aruguments:
        temp_a (ndarray) : Ambient temperature; degree Celsius
        speed (ndarray) : Wind speed; meters/second
        pres (ndarray) : Atmospheric pressure; hPa
        solar (ndarray) : Radiant heat flux incident on the globe;
            currently assuming this to be the solar irradiance; W/m**2
  
    Keyword arguments:
  
    Returns:
        ndarray : Globe temperature; degree Celsius
  
    """

    temp_g = numpy.empty( temp_a.shape, dtype=numpy.float32 )
    cdef:
        Py_ssize_t i, size = temp_a.size
        float [::1] temp_g_view = temp_g 
        float [::1] temp_a_view = temp_a.astype( numpy.float32 )
        float [::1] speed_view  = speed.astype(  numpy.float32 )
        float [::1] pres_view   = pres.astype(   numpy.float32 )
        float [::1] solar_view  = solar.astype(  numpy.float32 )
        float [::1] f_db_view   = f_db.astype(   numpy.float32 )
        float [::1] cosz_view   = cosz.astype(   numpy.float32 )


    for i in prange( size, nogil=True ):
        temp_g_view[i] = _globe_temperature( 
            temp_a_view[i]+CtoK,
            speed_view[i],
            pres_view[i],
            solar_view[i],
            f_db_view[i],
            cosz_view[i],
        ) - CtoK
    return temp_g

def psychrometric_wetbulb(temp_a, vapor_a=None, temp_d=None, relhum=None ):
    """
    Calculate psychrometric wet bulb from other variables

    Arguments:
        temp_a (ndarray) : ambient temperature in degree Celsius

    Keyword arguments:
        vapor_a (ndarray) : ambient vapor pressure in kiloPascals
        temp_d (ndarray) : Dew point temperature in degree Celsius
        relhum (ndarray) : Relative humidity as fraction

    """

    if vapor_a is None:
        if temp_d is not None:
            vapor_a = saturation_vapor_pressure( temp_d )/10.0
        elif relhum is not None:
            vapor_a = relhum * saturation_vapor_pressure( temp_a )/10.0
        else:
            raise Exception( 'Must input one of "vapor_a", "relhum", or "temp_d"' )

    return 0.376 + 5.79*vapor_a + (0.388 - 0.0465*vapor_a)*temp_a

def natural_wetbulb( temp_a, temp_psy, temp_g, speed ):
    """
    Compute natural wet bulb temperature using Bernard method
  
    If the globe temperature exceeds the dry bulb temperature by more than 
    4 degree C, then the effect of radiant heat is considered in the 
    calculation (factor e is computed). If this difference is 4 degree C or
    less, then no radiant heat effect is considered (factor C is computed)
  
    Arguments:
        temp_a (ndarray) : Ambient temperature; degree Celsius
        temp_psy (ndarray) : Psychrometric wet bulb temperature; degree Celsius
        temp_g (ndarray) : Globe temperature; degree Celsius
        speed (ndarray) : Wind speed; meters/second
  
    Keyword arguments:
        None.
  
    Returns:
        ndarray : Natural wet bulb temperature in degree Celsius
  
    Method obtained from:
        https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GH000231
  
    Notes:
        In the reference paper, the psychrometric wet bulb temperture is 
        Tpwb, however, to be consistent with other papers, and because I
        believe it makes more sense, we use temp_psy for psychrometric wet bulb
  
    """

    return numpy.where( (temp_g-temp_a) < 4.0,
      temp_a - factor_c( speed )*(temp_a - temp_psy),
      temp_psy + 0.25*(temp_g-temp_a) + factor_e( speed )
    )

def bernard( lat, lon, datetime,
        solar, pres, temp_a, temp_d, speed, 
        f_db   = None,
        cosz   = None,
        zspeed = units.Quantity( 10.0, 'meter' ),
        **kwargs):
    """
    Compute WBGT using Bernard Method

    Arguments:
        lat (float) : Latitude of observations
        lon (float) : Longitude of observations
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        solar (Quantity) : Solar irradiance; unit of power over area
        pres (Quantity) : Atmospheric pressure; unit of pressure
        temp_a (Quantity) : Ambient temperature; unit of temperature
        temp_d (Quantity) : Dew point temperature; unit of temperature
        speed (Quantity) : Wind speed; units of speed

    Keyword arguments:
        f_db (float) : Direct beam radiation from the sun; fraction
        cosz (float) : Cosine of solar zenith angle
        zspeed (Quantity) : Height of the wind speed measurment

    Returns:
        dict : 
          - Tg : Globe temperatures in ndarray
          - Tpsy : psychrometric wet bulb temperatures in ndarray
          - Tnwb : Natural wet bulb temperatures in ndarray
          - Twbg : Wet bulb-globe temperatures in ndarray

    """

    if (f_db is None) or (cosz is None):
        solar = solar_parameters( 
            lat, lon, datetime, solar.to('watt/m**2').magnitude, **kwargs,
        )
        if cosz is None:
            cosz = solar[1]
        if f_db is None:
            f_db = solar[2]
        solar = solar[0]

    temp_a = temp_a.to( 'degree_Celsius' ).magnitude
    temp_d = temp_d.to( 'degree_Celsius' ).magnitude
    pres   = pres.to(   'hPa'            ).magnitude
    speed  = loglaw( speed, zspeed, ).to('meter/second').magnitude

    temp_g = globe_temperature(
        temp_a, speed, pres, solar, f_db, cosz,
    )
    temp_psy = psychrometric_wetbulb( temp_a, temp_d = temp_d )
    temp_nwb = natural_wetbulb( temp_a, temp_psy, temp_g, speed )

    return {
        'Tg'    : units.Quantity(temp_g,   'degree_Celsius'),
        'Tpsy'  : units.Quantity(temp_psy, 'degree_Celsius'),
        'Tnwb'  : units.Quantity(temp_nwb, 'degree_Celsius'),
        'Twbg'  : units.Quantity(0.7*temp_nwb + 0.2*temp_g + 0.1*temp_a, 'degree_Celsius'),
        'solar' : units.Quantity( solar, 'watt/m**2' ),
    }
