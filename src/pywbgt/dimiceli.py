"""
WBGT from the Dimiceli method

Note that when running the main wetbulb_globe() function, wind
speed is adjusted to the 2 meter height and solar irradiance
may be adjusted based on the Liljegren method for calculating
solar parameters. Running other functions such as
globe_temperature() or natural_wetbulb calculations without
these adjustments may lead to different results.

"""

import numpy as np
from metpy.units import units

from .constants import SIGMA, MIN_SPEED, DIMICELI_MIN_SPEED
from .liljegren import solar_parameters
from .psychrometric_wetbulb import stull
from .natural_wetbulb import malchaire, hunter_minyard, nws_boyer
from .calc import relative_humidity, loglaw
from .utils import datetime_check

def adjust_speed_2m(speed, zspeed=None, min_speed=MIN_SPEED):
    """
    Adjust to 2m wind speed and clip

    Given the height of the wind speed measurment, adjust
    to the 2 meter wind speed. Also, esure that wind speeds
    are not less than min_speed.

    Arguments:
        speed (Quantity) : Wind speed; units of speed

    Keyword arguments:
        zspeed (Quantity) : Height of the wind speed measurment.
            Default is 10 meters
        min_speed (Quantity) : Sets the minimum speed for the height-adjusted
            wind speed. If this keyword is set, the larger of input value and
             DIMICELI_MIN_SPEED is used. The default value is MIN_SPEED, which
            is 2 knots, and is larger than DIMICELI_MIN_SPEED.
 

    Returns:
        tuple : Wind speed adjusted to 2 meter value and the
            minimum wind speed that the adjusted values have
            been clipped to.

    """

    if zspeed is None:
        zspeed = units.Quantity(10.0, 'meter')

    min_speed = max(min_speed, DIMICELI_MIN_SPEED)
    speed2m   = np.clip(
        loglaw( speed, zspeed ),
        min_speed,
        None,
    )

    return speed2m, min_speed

def conv_heat_flow_coeff(solar=None, zenith=None, **kwargs):
    """
    Convective heat flow coefficient

    Keyword arguments:
        solar (ndarray) : Solar irradiance
        zenith (ndarray) : Solar zenith angle (radians)

    Returns:
        float : The convective heat flow coefficient

    Notes:
        Default value obtained from:
        https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
  
    """

    const_a = const_b = const_c = 0.0
    if (solar is None) or (zenith is None):
        return 0.315

    return (
        const_a * 
        solar**const_b *
        np.cos(zenith)**const_c
    )


def atmospheric_vapor_pressure( temp_air, temp_dew, pres ):
    """
    Compute atmospheric vapor pressure

    Arguments:
        temp_air (float) : ambient temperature in degrees Celsius
        temp_dew (float) : dew point temperature in degrees Celsius
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
        (temp_d-temp_air) is odd. It is some kind of differential pressure factor.
        Some quick testing indicates that this formula will ALWAYS give
        a slighly lower pressure value than the Bolton (1980) formula; however,
        the values tend to be farily similar as long as there is not a large
        differenc between temp_air and temp_d

    """

    return (
        np.exp( (17.67 * (temp_dew - temp_air) ) / (temp_dew + 243.5) ) *
        (1.0007 + 3.46e-6 * pres) *
        6.112 * np.exp( 17.502 * temp_air / (240.97 + temp_air) )
    )

def thermal_emissivity( temp_air, temp_dew, pres ):
    """
    Compute thermal emissivity from readily available NWS values
  
    Arguments:
        temp_air (float) : ambient temperature in degrees Celsius
        temp_dew (float) : dew point temperature in degrees Celsius
        pres (float) : Barometric pressure in hPa
  
    Returns:
        float : thermal emissivity
  
    """

    return (
        0.575 *
        atmospheric_vapor_pressure( temp_air, temp_dew, pres )**(1.0/7.0)
    )

def factor_b( temp_air, temp_dew, pres, solar, f_db, cosz ):
    """

    Arguments:
        temp_air (float) : ambient temperature in degrees Celsius
        temp_dew (float) : dew point temperature in degrees Celsius
        pres (float) : Barometric pressure in hPa 
        solar (float) : solar irradiance in Watts per meter**2
        f_db (float) : Fraction of direct beam radiation
        cosz (float) : Cosine of solar zenith angle
  
    """

    f_dif = 1.0 - f_db
    return (
        solar * ( f_db/(4.0*SIGMA*cosz) + 1.2*f_dif/SIGMA ) +
        thermal_emissivity( temp_air, temp_dew, pres ) * temp_air**4
    )

def factor_c(speed, chfc=None, **kwargs):
    """
  
    Arguments:
        speed (float) : wind speed in meters per hour adjusted to 2 meters
            and clipped to (at the smallest value) DIMICELI_MIN_SPEED
  
    Keyword arguments:
        chfc (ndarray) : Convective heat flow coefficient. If none
            provided, is computed using the conv_heat_flow_coeff()
            function
        **kwargs : Any extra arguments are passed directly to
            the conv_heat_flow_coeff() function.

    """

    if chfc is None:
        chfc = conv_heat_flow_coeff(**kwargs)
    return chfc * speed**0.58 / 5.3865e-8

def globe_temperature( temp_air, temp_dew, pres, speed, solar, f_db, cosz ):
    """
    Compute globe temperature
  
    Arguments:
        temp_air (float) : ambient temperature in degrees Celsius
        temp_dew (float) : dew point temperature in degrees Celsius
        pres (float) : Barometric pressure in hPa 
        speed (float) : wind speed in meters per hour adjusted to
            2 meters and clipped to (at the smallest value) 
            DIMICELI_MIN_SPEED
        solar (float) : solar irradiance in Watts per meter**2
        f_db (float) : Fraction of direct beam radiation
        cosz (float) : Cosine of solar zenith angle
 
    Returns:
        ndarray : Black globe temperature in degrees C

    Notes:
        Chapter 26 of IAENG Transactions on Engineering Technologies: 
          "Black Globe Temperature Estimate for the WBGT Index"
    
        https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
  
    """

    fac_b = factor_b( temp_air, temp_dew, pres, solar, f_db, cosz)
    fac_c = factor_c( speed )

    return (fac_b + fac_c*temp_air + 7.68e6) / (fac_c + 2.56e5)

def psychrometric_wetbulb( temp_air, temp_dew ):
    """
    Wet bulb temperature from Dimiceli method
  
    This formula for wet bulb temperature appears at the bottom of 
    "Estimation of Black Globe Temperature for Calculation of the WBGT Index"
    by Dimiceli and Piltz.
  
    https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
  
    Inputs:
        temp_air (ndarray) : Ambient (dry bulb) temperature (degree C)
        temp_dew (ndarray) : Dew point temperature (degree C)
  
    """

    relhum = relative_humidity( temp_air, temp_dew ) * 100.0
    return (
           -5.806    + 0.672   *temp_air -  0.006   *temp_air**2       +
         (  0.061    + 0.004   *temp_air + 99.000e-6*temp_air**2) * relhum +
         (-33.000e-6 - 5.000e-6*temp_air -  1.000e-7*temp_air**2) * relhum**2
    )

def wetbulb_globe(
        datetime, lat, lon,
        solar, pres, temp_air, temp_dew, speed,
        f_db      = None,
        cosz      = None,
        zspeed    = None,
        min_speed = MIN_SPEED,
        **kwargs,
    ):
    """
    Compute WBGT using Dimiceli method

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
             DIMICELI_MIN_SPEED is used. The default value is MIN_SPEED, which
            is 2 knots, and is larger than DIMICELI_MIN_SPEED.
        wetbulb (str) : Name of wet bulb algorithm to use:
            {dimiceli, stull} DEFAULT = dimiceli
        natural_wetbulb (str) : Name of the natural wetbulb algorithm to use:
            (hunter_minyard, malchaire, boyer). Default is hunter_minyard.

    Notes:
        Enacts 'Rectent Updates and Improvements' from the following white paper:

            https://vlab.noaa.gov/documents/6609493/7858379/NDFD+WBGT+Description+Document.pdf/fb89cc3a-0536-111f-f124-e4c93c746ef7?t=1642792547129

    Returns:
        dict : 
            - Tg : Globe temperatures in ndarray
            - Tpsy : psychrometric wet bulb temperatures in ndarray
            - Tnwb : Natural wet bulb temperatures in ndarray
            - Twbg : Wet bulb-globe temperatures in ndarray
            - solar : Solar irradiance from Liljegren

    """

    datetime = datetime_check(datetime)

    solar     = solar.to(   'watt/m**2'     ).magnitude
    pres      = pres.to(    'hPa'           ).magnitude
    temp_air  = temp_air.to('degree_Celsius').magnitude
    temp_dew  = temp_dew.to('degree_Celsius').magnitude
    speed2m, min_speed = adjust_speed_2m(
        speed,
        zspeed    = zspeed,
        min_speed = min_speed,
    )

    if (f_db is None) or (cosz is None):
        solar = solar_parameters(datetime, lat, lon, solar, **kwargs )
        if cosz is None:
            cosz = solar[1]
        if f_db is None:
            f_db = solar[2]
        solar = solar[0]

    temp_g = globe_temperature(
        temp_air,
        temp_dew,
        pres,
        speed2m.to('meter per hour').magnitude,
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

    nwb_method = kwargs.get('natural_wetbulb', 'HUNTER_MINYARD').upper()
    if nwb_method == 'HUNTER_MINYARD':
        temp_nwb  = hunter_minyard(
            temp_psy,
            solar*f_db,
            speed2m.to('meter per second').magnitude,
        )
    elif nwb_method == 'MALCHAIRE':
        temp_nwb  = malchaire( temp_air, temp_dew, temp_psy, temp_g )
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
        'Tg'        : units.Quantity(temp_g, 'degree_Celsius'),
        'Tpsy'      : units.Quantity(temp_psy, 'degree_Celsius'),
        'Tnwb'      : units.Quantity(temp_nwb, 'degree_Celsius'), 
        'Twbg'      : units.Quantity(0.7*temp_nwb + 0.2*temp_g + 0.1*temp_air, 'degree_Celsius'),
        'solar'     : units.Quantity( solar, 'watt/m**2'),
        'speed'     : speed2m.to('meter per second'),
        'min_speed' : min_speed,
    }
