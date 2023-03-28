"""
WBGT from the Dimiceli method

"""

import numpy
from metpy.units import units

from . import SIGMA
from .calc import relative_humidity, loglaw
from .liljegren import solar_parameters
from .psychrometric_wetbulb import stull
from .natural_wetbulb import malchaire, hunter_minyard, nws_boyer

def conv_heat_flow_coeff( solar=None, zenith=None):
    """
    Convective heat flow coefficient

    Keyword arguments:
        solar () : Solar irradiance
        zenith () : Solar zenith angle (radians)

    Returns:
        float : The convective heat flow coefficient

    Notes:
        Default value obtained from:
        https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
  
    """

    a = b = c = 0.0
    if (solar is not None) and (zenith is not None):
        return a * solar**b * numpy.cos( zenith )**c

    return 0.315

def atmospheric_vapor_pressure( temp_a, temp_d, pres ):
    """
    Compute atmospheric vapor pressure

    Arguments:
        temp_a (float) : ambient temperature in degrees Celsius
        temp_d (float) : dew point temperature in degrees Celsius
        pres (float) : Barometric pressure in hPa

    Returns:
        float : atmospheric vapor pressure

    Note:
        This is an odd function that I cannot find a lot of informaiton
        about. I was able to determe that the 2nd and 3rd lines of the function
        come from Buck, A. 1981: New Equations for computing vapor pressure
        and enhancement factor and is an improved equation for calculating
        vapor pressure in air, as opposed to pure water vapor.
  
        However, the first line of the equation does not make much sense.
        The base of the equation is from Bolton (1980), but the use of
        (temp_d-temp_a) is odd. It is some kind of differential pressure factor.
        Some quick testing indicates that this formula will ALWAYS give
        a slighly lower pressure value than the Bolton (1980) formula; however,
        the values tend to be farily similar as long as there is not a large
        differenc between temp_a and temp_d

    """

    return numpy.exp( (17.67 * (temp_d - temp_a) ) / (temp_d + 243.5) ) *\
      (1.0007 + 3.46e-6 * pres) *\
      6.112 * numpy.exp( 17.502 * temp_a / (240.97 + temp_a) )

def thermal_emissivity( temp_a, temp_d, pres ):
    """
    Compute thermal emissivity from readily available NWS values
  
    Arguments:
        temp_a (float) : ambient temperature in degrees Celsius
        temp_d (float) : dew point temperature in degrees Celsius
        pres (float) : Barometric pressure in hPa
  
    Returns:
        float : thermal emissivity
  
    """

    return 0.575 * atmospheric_vapor_pressure( temp_a, temp_d, pres )**(1.0/7.0)

def factor_b( temp_a, temp_d, pres, solar, f_db, cosz ):
    """

    Arguments:
        temp_a (float) : ambient temperature in degrees Celsius
        temp_d (float) : dew point temperature in degrees Celsius
        pres (float) : Barometric pressure in hPa 
        solar (float) : solar irradiance in Watts per meter**2
        f_db (float) : Fraction of direct beam radiation
        cosz (float) : Cosine of solar zenith angle
  
    """

    f_dif = 1.0 - f_db
    return solar * ( f_db/(4.0*SIGMA*cosz) + 1.2*f_dif/SIGMA ) +\
        thermal_emissivity( temp_a, temp_d, pres ) * temp_a**4

def factor_c( speed, chfc = None ):
    """
  
    Arguments:
      speed (float) : wind speed in meters per hour
  
    """

    if chfc is None:
        chfc = conv_heat_flow_coeff()
    return chfc * speed**0.58 / 5.3865e-8

def globe_temperature( temp_a, temp_d, pres, speed, solar, f_db, cosz ):
    """
    Compute globe temperature
  
    Arguments:
        temp_a (float) : ambient temperature in degrees Celsius
        temp_d (float) : dew point temperature in degrees Celsius
        pres (float) : Barometric pressure in hPa 
        speed (float) : wind speed in meters per hour
        solar (float) : solar irradiance in Watts per meter**2
        f_db (float) : Fraction of direct beam radiation
        cosz (float) : Cosine of solar zenith angle
  
    Notes:
        Chapter 26 of IAENG Transactions on Engineering Technologies: 
          "Black Globe Temperature Estimate for the WBGT Index"
    
        https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
  
    """

    fac_b = factor_b( temp_a, temp_d, pres, solar, f_db, cosz)
    fac_c = factor_c( speed )

    return (fac_b + fac_c*temp_a + 7.68e6) / (fac_c + 2.56e5)

def psychrometric_wetbulb( temp_a, temp_d ):
    """
    Wet bulb temperature from Dimiceli method
  
    This formula for wet bulb temperature appears at the bottom of 
    "Estimation of Black Globe Temperature for Calculation of the WBGT Index"
    by Dimiceli and Piltz.
  
    https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
  
    Inputs:
        temp_a (ndarray) : Ambient (dry bulb) temperature (degree C)
        temp_d (ndarray) : Dew point temperature (degree C)
  
    """

    relhum = relative_humidity( temp_a, temp_d ) * 100.0
    return -5.806    + 0.672   *temp_a -  0.006   *temp_a**2       +\
         (  0.061    + 0.004   *temp_a + 99.000e-6*temp_a**2) * relhum +\
         (-33.000e-6 - 5.000e-6*temp_a -  1.000e-7*temp_a**2) * relhum**2

def dimiceli( lat, lon, datetime,
              solar, pres, temp_air, temp_dew, speed,
              f_db   = None,
              cosz   = None,
              zspeed = units.Quantity(10.0, 'meter'),
              **kwargs ):
    """
    Compute WBGT using Dimiceli method

    Arguments:
        lat (float) : Latitude of observations
        lon (float) : Longitude of observations
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        solar (Quantity) : Solar irradiance; unit of power over area
        pres (Quantity) : Atmospheric pressure; unit of pressure
        temp_air (Quantity) : Ambient temperature; unit of temperature
        temp_dew (Quantity) : Dew point temperature; unit of temperature
        speed (Quantity) : Wind speed; units of speed

    Keyword arguments:
        f_db (float) : Direct beam radiation from the sun; fraction
        cosz (float) : Cosine of solar zenith angle
        zspeed (Quantity) : Height of the wind speed measurment
        wetbulb (str) : Name of wet bulb algorithm to use:
          {dimiceli, stull} DEFAULT = dimiceli
        natural_wetbulb (str) : Name of the natural wet bulb algorithm to use:
          {malchaire, hunter_minyard} DEFAULT = malchaire

    Returns:
        dict : 
          - Tg : Globe temperatures in ndarray
          - Tpsy : psychrometric wet bulb temperatures in ndarray
          - Tnwb : Natural wet bulb temperatures in ndarray
          - Twbg : Wet bulb-globe temperatures in ndarray
  
    """

    solar    = solar.to( 'watt/m**2'       ).magnitude
    pres     = pres.to( 'hPa'             ).magnitude
    temp_air = temp_air.to( 'degree_Celsius'  ).magnitude
    temp_dew = temp_dew.to( 'degree_Celsius'  ).magnitude
    speed2m  = loglaw( speed, zspeed )
    # Ensure wind speed is at least one (1) mile/hour
    speed1   = numpy.clip(speed2m.to('meters per hour').magnitude, 1690.0, None)

    if (f_db is None) or (cosz is None):
        solar = solar_parameters( lat, lon, datetime, solar, **kwargs )
        if cosz is None:
            cosz = solar[1]
        if f_db is None:
            f_db = solar[2]
        solar = solar[0]

    temp_g = globe_temperature(
        temp_air,
        temp_dew,
        pres,
        speed1,
        solar,
        f_db,
        cosz,
    )

    wb_method = kwargs.get('wetbulb', 'DIMICELI').upper()
    if wb_method == 'DIMICELI':
        temp_psy = psychrometric_wetbulb( temp_air, temp_dew)
    elif wb_method == 'STULL':
        temp_psy = stull( temp_air, temp_dew)
    else:
        raise Exception( f"Invalid option for 'wetbulb' : {wb_method}" )

    nwb_method = kwargs.get('natural_wetbulb', 'MALCHAIRE').upper()
    if nwb_method == 'MALCHAIRE':
        temp_nwb  = malchaire( temp_air, temp_dew, temp_psy, temp_g )
    elif nwb_method == 'HUNTER_MINYARD':
        temp_nwb  = hunter_minyard(
            temp_psy,
            solar*f_db,
            speed2m.to('meter per second').magnitude,
        )
    elif nwb_method == 'BOYER':
        temp_nwb  = nws_boyer(
            temp_air,
            temp_psy,
            solar*f_db,
            speed2m.to('meter per second').magnitude,
        )
    else:
        raise Exception(
            f"Unsupported value for 'natural_wetbulb' : {nwb_method}"
        )

    return {
        'Tg'    : units.Quantity(temp_g, 'degree_Celsius'),
        'Tpsy'  : units.Quantity(temp_psy, 'degree_Celsius'),
        'Tnwb'  : units.Quantity(temp_nwb, 'degree_Celsius'), 
        'Twbg'  : units.Quantity(0.7*temp_nwb + 0.2*temp_g + 0.1*temp_air, 'degree_Celsius'),
        'solar' : units.Quantity( solar, 'watt/m**2'),
    }
