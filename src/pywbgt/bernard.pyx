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

from libc.math cimport fabs, pow, log10
cimport numpy
cimport cython

from metpy.units import units

from .constants import SIGMA, EPSILON, MIN_SPEED
from .liljegren import solar_parameters
from .calc import saturation_vapor_pressure, loglaw
from .utils import datetime_check

cdef:
    int   MAX_ITER  = 50
    float CtoK      = 273.15
    float CONVERGE  = 1.0E-3
    float ALPHA_SFC = 0.00
    float NaN       = numpy.nan
    float SIGMAB    = SIGMA
    float EPS       = EPSILON

@cython.cdivision(True)
cdef double _conv_heat_trans_coeff(
        double temp_g, double temp_air, double speed,
    ) nogil:
    """
    Convective heat transfer coefficient

    Compute the convective heat transfer coefficient
    from Bernard method

    Arguments:
        temp_g : Globe temperature; Kelvin
        temp_air : Air temperature; Kelvin
        speed : Wind speed; meter/second

    Returns:
        convective heat transfer coefficient

    """

    cdef double coeff, delta_t = temp_g-temp_air
    coeff = pow(
        pow(10.9*pow(speed, 0.566), 3) + 
        pow(0.35 + 1.77*pow(fabs(delta_t), 0.25), 3),
        1.0/3.0
    )
    if (delta_t < 0.0): 
        return -coeff
    return coeff

@cython.cdivision(True)
cdef double _globe_temperature(
        double temp_air,
        double speed,
        double pres,
        float solar, 
        float f_db,
        float cosz,
    ) nogil:
    """
    Determine globe temperature through iterative solver

    Aruguments:
        temp_air (ndarray) : Ambient temperature; Kelvin 
        speed (ndarray) : Wind speed; meters/second
        pres (ndarray) : Atmospheric pressure; hPa
        solar (ndarray) : Radiant heat flux incident on the globe;
            currently assuming this to be the solar irradiance; W/m**2

    Keyword arguments:
        None.

    Returns:
        ndarray : Globe temperature determined through iterative solver

    """

    cdef:
        int ii
        double h, temp_g_new, temp_g = temp_air

    for ii in range( MAX_ITER ):
        h = _conv_heat_trans_coeff( temp_g, temp_air, speed )
        temp_g_new = pow(
            temp_air**4 +
            solar/SIGMAB/2.0 * (1 + f_db*(1/2.0/cosz - 1.0) + ALPHA_SFC) -
            h/EPS/SIGMAB * (temp_g-temp_air),
            0.25
        )

        if fabs(temp_g_new-temp_g) < CONVERGE: 
            return temp_g_new
        temp_g = 0.9*temp_g + 0.1*temp_g_new

    return NaN

def conv_heat_trans_coeff( temp_g, temp_air, speed ):
    """
    Convective heat transfer coefficent
  
    This is equation 5 from Bernard and Pourmoghani (1999) used to
    compute the convective heat transfer coefficient.
  
    Arguments:
        temp_g (ndarray) : Globe temperature; degree Celsius
        temp_air (ndarray) : Ambient temperature; degree Celsius
        speed (ndarray) : Wind speed; meters/second
  
    Kewyord arguments:
        None.
  
    Returns:
        ndarray : Convective heat transfer coefficient values
  
    """

    delta_t = temp_g-temp_air
    coeff = (
        (10.9*speed**0.566)**3 + (0.35 + 1.77*abs(delta_t)**0.25)**3
    )**(1.0/3.0)

    return numpy.where(
        delta_t < 0,
        -coeff,
        coeff,
    )

@cython.cdivision(True)
cdef double _factor_c(double speed) nogil:
    """
    Compute the C factor on a single value

    Arguments:
        speed (double) : Wind speed in meters per second

    Returns:
        double : Value of the C factor.

    """

    if speed < 0.03:
        return 0.85
    if speed > 3.0:
        return 1.0
    return 0.96 + 0.069*log10(speed)
  
def factor_c( speed ):
    """
    Compute factor for natural wet bulb temperature without radiant heat
  
    This factor is used to relate psychrometric wet bulb temperature
    to natural wet bulb temperature in the absence of radiant heat.
  
    Arguments:
        speed (ndarray) : Wind speed in meters/second
  
    """

    fac_c      = numpy.full( speed.shape, 0.85 )
    idx        = numpy.where( speed>= 0.03 )
    fac_c[idx] = 0.96 + 0.069*numpy.log10( speed[idx] )
    # Where wind > 3.0, keep values of C, else compute C and return values
    return numpy.where( speed > 3.0, 1.0, fac_c )

@cython.cdivision(True)
cdef double _factor_e(double speed) nogil:
    """
    Compute the E factor on a single value

    Arguments:
        speed (double) : Wind speed in meters per second

    Returns:
        double : Value of the E factor.

    """

    if speed < 0.1:
        return 1.1
    if speed > 1.0:
        return -0.1
    return 0.1/pow(speed, 1.1) - 0.2
 
def factor_e( speed ):
    """
    Compute factor for natural wet bulb temperature with radiant heat
  
    This factor is used to related psychrometric wet bulb temperature
    to natural wet bulb temperature in the presence of radiant heat.
  
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
def globe_temperature_64(
        double [::1] temp_air,
        double [::1] speed,
        double [::1] pres,
        float  [::1] solar,
        float  [::1] f_db,
        float  [::1] cosz,
    ): 
    """
    Compute Tg (64-bit)

    Compute value(s) for globe temperature in double precision

    """

    temp_g = numpy.empty(temp_air.size, dtype=numpy.float64)
    cdef:
        Py_ssize_t i, size = temp_air.size
        double [::1] temp_g_view   = temp_g 

    for i in prange( size, nogil=True ):
        temp_g_view[i] = _globe_temperature( 
            temp_air[i]+CtoK,
            speed[i],
            pres[i],
            solar[i],
            f_db[i],
            cosz[i],
        ) - CtoK
    return temp_g

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def globe_temperature_32(
        float [::1] temp_air,
        float [::1] speed,
        float [::1] pres,
        float [::1] solar,
        float [::1] f_db,
        float [::1] cosz,
    ): 
    """
    Compute Tg (32-bit)

    Compute value(s) for globe temperature in single precision

    """

       
    temp_g = numpy.empty(temp_air.size, dtype=numpy.float32)
    cdef:
        Py_ssize_t i, size = temp_air.size
        float [::1] temp_g_view = temp_g 

    for i in prange( size, nogil=True ):
        temp_g_view[i] = <float>_globe_temperature( 
            <double>temp_air[i]+CtoK,
            <double>speed[i],
            <double>pres[i],
            solar[i],
            f_db[i],
            cosz[i],
        ) - CtoK
    return temp_g

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def globe_temperature(temp_air, speed, pres, solar, f_db, cosz):
    """
    Determine globe temperature through iterative solver

    Arguments:
        temp_air (ndarray) : Ambient temperature; degree Celsius
        speed (ndarray) : Wind speed; meters/second
        pres (ndarray) : Atmospheric pressure; hPa
        solar (ndarray) : Radiant heat flux incident on the globe;
            currently assuming this to be the solar irradiance; W/m**2

    Returns:
        ndarray : Globe temperature; degree Celsius

    """

    # if these variables are NOT all the same type, make them all float32
    if not temp_air.dtype == speed.dtype == pres.dtype:
        temp_air = temp_air.astype(numpy.float32)
        speed    =    speed.astype(numpy.float32)
        pres     =     pres.astype(numpy.float32)

    # If these variables are NOT all float32, force to float32
    if not solar.dtype == f_db.dtype == cosz.dtype == numpy.float32:
        solar = solar.astype(numpy.float32)
        f_db  =  f_db.astype(numpy.float32)
        cosz  =  cosz.astype(numpy.float32)

    # Run 64-bit version
    if temp_air.dtype == numpy.float64:
        return globe_temperature_64(temp_air, speed, pres, solar, f_db, cosz)
    # Run 32-bit version
    if temp_air.dtype == numpy.float32:
        return globe_temperature_32(temp_air, speed, pres, solar, f_db, cosz)
    # Error as MUST input floating-point values
    raise Exception('Must imput floating-point values')


def psychrometric_wetbulb(temp_air, vapor_air=None, temp_dew=None, relhum=None):
    """
    Calculate psychrometric wet bulb from other variables

    Arguments:
        temp_air (ndarray) : ambient temperature in degree Celsius

    Keyword arguments:
        vapor_air (ndarray) : ambient vapor pressure in kiloPascals
        temp_dew (ndarray) : Dew point temperature in degree Celsius
        relhum (ndarray) : Relative humidity as fraction

    """

    if vapor_air is None:
        if temp_dew is not None:
            vapor_air = saturation_vapor_pressure( temp_dew )/10.0
        elif relhum is not None:
            vapor_air = relhum * saturation_vapor_pressure( temp_air )/10.0
        else:
            raise Exception(
                'Must input one of "vapor_air", "relhum", or "temp_dew"'
            )

    return 0.376 + 5.79*vapor_air + (0.388 - 0.0465*vapor_air)*temp_air

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def natural_wetbulb_64(
        double [::1] temp_air,
        double [::1] temp_psy,
        double [::1] temp_g,
        double [::1] speed,
    ):
    """
    Compute Tnwb (64-bit)

    Compute value(s) for natural wetbulb temperature in double precision

    """
 
    temp_nwb = numpy.empty(temp_air.size, dtype=numpy.float64)
    cdef:
        Py_ssize_t i, size = temp_air.size
        double val 
        double [::1] temp_nwb_view = temp_nwb

    for i in prange( size, nogil=True ):
        val = temp_g[i]-temp_air[i]
        if val < 4.0:
            val = temp_air[i] - temp_psy[i]
            val = temp_air[i] - _factor_c(speed[i]) * val
        else:
            val = temp_psy[i] + 0.25*val + _factor_e(speed[i])

        temp_nwb_view[i] = val
    return temp_nwb

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def natural_wetbulb_32(
        float [::1] temp_air,
        float [::1] temp_psy,
        float [::1] temp_g,
        float [::1] speed,
    ):
    """
    Compute Tnwb (32-bit)

    Compute value(s) for natural wetbulb temperature in single precision

    """
  
    temp_nwb = numpy.empty(temp_air.size, dtype=numpy.float32)
    cdef:
        Py_ssize_t i, size = temp_air.size
        float val 
        float [::1] temp_nwb_view = temp_nwb

    for i in prange( size, nogil=True ):
        val = temp_g[i]-temp_air[i]
        if val < 4.0:
            val = temp_air[i] - temp_psy[i]
            val = temp_air[i] - <float>_factor_c(<double>speed[i]) * val
        else:
            val = temp_psy[i] + 0.25*val + <float>_factor_e(<double>speed[i])

        temp_nwb_view[i] = val
    return temp_nwb

def natural_wetbulb( temp_air, temp_psy, temp_g, speed ):
    """
    Compute natural wet bulb temperature using Bernard method
  
    If the globe temperature exceeds the dry bulb temperature by more than 
    4 degree C, then the effect of radiant heat is considered in the 
    calculation (factor e is computed). If this difference is 4 degree C or
    less, then no radiant heat effect is considered (factor C is computed)
  
    Arguments:
        temp_air (ndarray) : Ambient temperature; degree Celsius
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

    if not temp_air.dtype == temp_psy.dtype == temp_g.dtype == speed.dtype:
        temp_air = temp_air.astype(numpy.float32)
        temp_psy = temp_psy.astype(numpy.float32)
        temp_g   =   temp_g.astype(numpy.float32)
        speed    =    speed.astype(numpy.float32)

    if temp_air.dtype == numpy.float64:
        return natural_wetbulb_64(temp_air, temp_psy, temp_g, speed)
    if temp_air.dtype == numpy.float32:
        return natural_wetbulb_32(temp_air, temp_psy, temp_g, speed)
    raise Exception('Must imput floating-point values')

def wetbulb_globe(
        datetime, lat, lon,
        solar, pres, temp_air, temp_dew, speed,
        f_db      = None,
        cosz      = None,
        zspeed    = None,
        min_speed = None, 
        **kwargs,
    ):
    """
    Compute WBGT using Bernard Method

    Arguments:
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        lat (float) : Latitude of observations
        lon (float) : Longitude of observations
        solar (Quantity) : Solar irradiance; unit of power over area
        pres (Quantity) : Atmospheric pressure; unit of pressure
        temp_air (Quantity) : Ambient temperature; unit of temperature
        temp_dew (Quantity) : Dew point temperature; unit of temperature
        speed (Quantity) : Wind speed; units of speed

    Keyword arguments:
        f_db (float) : Direct beam radiation from the sun; fraction
        cosz (float) : Cosine of solar zenith angle
        zspeed (Quantity) : Height of the wind speed measurment.
            Default is 10 meters
        min_speed (Quantity) : Sets the minimum speed for the height-adjusted
            wind speed. If this keyword is set, the larger of input value and
            MIN_SPEED is used. The default value is MIN_SPEED, which
            is 2 knots.

    Returns:
        dict : 
            - Tg : Globe temperatures as Quantity
            - Tpsy : psychrometric wet bulb temperatures as Quantity
            - Tnwb : Natural wet bulb temperatures as Quantity
            - Twbg : Wet bulb-globe temperatures as Quantity
            - solar : Solar irradiance from Liljegren as Quantity
            - speed : 2 meter adjusted wind speed as Quantity; will be same as input if already 2m wind speed
            - min_speed : Minimum speed that adjusted wind speed is clipped to as Quantity

    """

    datetime = datetime_check(datetime)
    if zspeed is None:
        zspeed = units.Quantity( 10.0, 'meter' )

    if (f_db is None) or (cosz is None):
        solar = solar_parameters( 
            datetime, lat, lon, solar.to('watt/m**2').magnitude, **kwargs,
        )
        if cosz is None:
            cosz = solar[1]
        if f_db is None:
            f_db = solar[2]
        solar = solar[0]

    temp_air = temp_air.to( 'degree_Celsius' ).magnitude
    temp_dew = temp_dew.to( 'degree_Celsius' ).magnitude
    pres   = pres.to(   'hPa'            ).magnitude

    if min_speed is None:
        min_speed = MIN_SPEED
    else:
        min_speed = max(
            min_speed,
            MIN_SPEED,
        )
 
    speed = numpy.clip(
        loglaw(speed, zspeed),
        min_speed,
        None,
    ).to('meter/second')

    temp_g = globe_temperature(
        temp_air, 
        speed.magnitude,
        pres,
        solar,
        f_db,
        cosz,
    )
    temp_psy = psychrometric_wetbulb(
        temp_air,
        temp_dew = temp_dew
    )
    temp_nwb = natural_wetbulb(
        temp_air,
        temp_psy,
        temp_g,
        speed.magnitude
    )

    return {
        'Tg'        : units.Quantity(temp_g,   'degree_Celsius'),
        'Tpsy'      : units.Quantity(temp_psy, 'degree_Celsius'),
        'Tnwb'      : units.Quantity(temp_nwb, 'degree_Celsius'),
        'Twbg'      : units.Quantity(0.7*temp_nwb + 0.2*temp_g + 0.1*temp_air, 'degree_Celsius'),
        'solar'     : units.Quantity( solar, 'watt/m**2' ),
        'speed'     : speed.to('meter/second'),
        'min_speed' : min_speed.to('meter/second'),
    }
