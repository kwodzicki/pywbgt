"""
Solar position and direct beam fraction

Python implementation of the Liljegren solar position and direct beam
fraction calculations utilizing the Python pvlib implementation of the
National Renewable Energy Laboratory (NREL) Solar Position Algorithm (SPA).

"""

import os
os.environ['PVLIB_USE_NUMBA'] = '1' # Force use of numba in pvlib.spa module

import numpy as np
from numba import njit, prange
from pvlib import spa

from .utils import datetime_adjust
from .liljegren import (
    LILJEGREN_CZA_MIN, LILJEGREN_NORMSOLAR_MAX, LILJEGREN_SOLAR_CONST,
)

# Default values for some parameters
ELEV     =    0.00
PRESSURE = 1013.25
TEMP     =   15.00

def solar_parameters(
    datetime, lat, lon, solar,
    gmt      = None,
    avg      = None,
    elev     = None,
    pressure = None,
    temp     = None,
     **kwargs
):
    """
    Calculate solar parameters based on date and location

    This function is adapted from the Liljegren WBGT C-code version 1.1 that
    is included in this package. As their solar position algorithm is only
    'valid' from 1950-2050, and we are unable to redistribute the NREL SPA
    source code, the Python pvlib package's implementation of the SPA algorithm
    is used.

    This function aims to replicated the funcationality of the 
    calc_solar_parameters() function from the Liljegren C-cdoe.

    The following quantities are calculated:
      - Modifies solar values to be consistent with normsolar from C function
      - cosine of the solar zenith angle
      - fraction of the solar irradiance due to the direct beam

    Arguments:
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        lat (ndarray) : Latitude of location(s) to compute parameters for; decimal
        lon (ndarray) : Longitude of location(s) to compute parameters for; decimal
        solar (ndarray) : Solar irradiance values (Watt/m**2)

    Keyword arguments:
        gmt (ndarray) : LST-GMT difference  (hours; negative in USA)
        avg (ndarray) : averaging time of the meteorological inputs (minutes)
        elev (ndarray) : Elevation of location (meters).
            Default is 0 m, or sea level.
        pressure (ndarray) : Average yearly pressure at location (hPa).
            Default is 1013.25 hPa
        temp (ndarray) : Average yearly temperature at location (degC)
            Default is 15.0 C

    Returns:
        tuple : Three (3) ndarrays containing:
            - Potentially modified solar radiation values
            - cosine of zenith angle
            - fraction of solar irradiance due to the direct beam
 
    """

    datetime = (
        datetime_adjust(datetime, gmt, avg)
        .values        
        .astype(np.int64)
    )/1.0e9

    ntime = datetime.shape[0]

    if lat.size <= 1:                                                         # If input latitude is only one (1) element, assume lon and urban are also one (1) element and expand all to match size of data
        lat = lat.repeat(ntime)

    if lon.size <= 1:
        lon = lon.repeat(ntime)

    if elev is None:
        elev = np.full(ntime, ELEV)
    elif elev.size <= 1:
        elev = elev.repeat(ntime)

    if pressure is None:
        pressure = np.full(ntime, PRESSURE)
    elif pressure.size <= 1:
        pressure = pressure.repeat(ntime)

    if temp is None:
        temp = np.full(ntime, TEMP)
    elif temp.size <= 1:
        temp = temp.repeat(ntime)
    
    solardist, zenith = solar_position(
        datetime,
        lat,
        lon,
        elev,
        pressure,
        temp,
        0,
        0,
    )

    cza      = np.cos(np.deg2rad(zenith))
    fdir     = np.zeros(cza.shape, dtype=np.float32)

    toasolar = LILJEGREN_SOLAR_CONST * np.clip(cza, 0.0, None) / solardist**2
    toasolar[cza < LILJEGREN_CZA_MIN] = 0.0

    idx      = toasolar > 0.0
    
    normsolar  = np.clip(solar[idx]/toasolar[idx], None, LILJEGREN_NORMSOLAR_MAX)
    solar[idx] = normsolar * toasolar[idx]
    fdir[idx]  = np.where(
        normsolar > 0.0,
        np.clip(np.exp(3.0-1.34*normsolar-1.65/normsolar), 0.0, 0.9),
        0.0,
    )

    return np.where(idx, solar, 0.0), cza, fdir

@njit(parallel=True)
def solar_position(
    unixtime,
    lat,
    lon,
    elev,
    pressure,
    temp,
    delta_t,
    atmos_refract,
): 
    """
    Adapted from pvlib.spa.solar_position_numpy

    This function is adapted from the numpy version of the numpy version]
    of the SPA algorithm provided by the pvlib package. This function has
    been wrapped as a numba parallel compute function to speed up computation.
    The returns from the function have also been modified as we are only
    interested in the Earth-Sun distance and the solar zenith angle.
    Some function calls from the solar_position_numpy() function have been
    removed to speed up computation as they are no longer requried.

    Arguments:
        See arguments for solar_position()

    Returns:
        tuple : numpy.ndarray for Earth-Sun distance (AU) and 
            solar zenith angle (degrees)

    """

    solardist = np.empty(unixtime.size, dtype=np.float32)
    zenith    = np.empty(unixtime.size, dtype=np.float32)

    for i in prange(unixtime.size):
        jd    = spa.julian_day(unixtime[i])
        jde   = spa.julian_ephemeris_day(jd, delta_t)
        jc    = spa.julian_century(jd)
        jce   = spa.julian_ephemeris_century(jde)
        jme   = spa.julian_ephemeris_millennium(jce)
        R     = spa.heliocentric_radius_vector(jme)
        solardist[i] = R

        L     = spa.heliocentric_longitude(jme)
        B     = spa.heliocentric_latitude(jme)
        Theta = spa.geocentric_longitude(L)
        beta  = spa.geocentric_latitude(B)
        x0    = spa.mean_elongation(jce)
        x1    = spa.mean_anomaly_sun(jce)
        x2    = spa.mean_anomaly_moon(jce)
        x3    = spa.moon_argument_latitude(jce)
        x4    = spa.moon_ascending_longitude(jce)
        
        l_o_nutation = np.empty((2,))
        spa.longitude_obliquity_nutation(jce, x0, x1, x2, x3, x4, l_o_nutation)

        delta_psi     = l_o_nutation[0]
        delta_epsilon = l_o_nutation[1]
        epsilon0      = spa.mean_ecliptic_obliquity(jme)
        epsilon       = spa.true_ecliptic_obliquity(epsilon0, delta_epsilon)
        delta_tau     = spa.aberration_correction(R)
        lamd          = spa.apparent_sun_longitude(Theta, delta_psi, delta_tau)
        v0            = spa.mean_sidereal_time(jd, jc)
        v             = spa.apparent_sidereal_time(v0, delta_psi, epsilon)
        alpha         = spa.geocentric_sun_right_ascension(lamd, epsilon, beta)
        delta         = spa.geocentric_sun_declination(lamd, epsilon, beta)

        H           = spa.local_hour_angle(v, lon[i], alpha)
        xi          = spa.equatorial_horizontal_parallax(R)
        u           = spa.uterm(lat[i])
        x           = spa.xterm(u, lat[i], elev[i])
        y           = spa.yterm(u, lat[i], elev[i])
        delta_alpha = spa.parallax_sun_right_ascension(x, xi, H, delta)
        delta_prime = spa.topocentric_sun_declination(
            delta, x, y, xi, delta_alpha, H,
        )
        H_prime     = spa.topocentric_local_hour_angle(H, delta_alpha)
        e0          = spa.topocentric_elevation_angle_without_atmosphere(
            lat[i], delta_prime, H_prime,
        )
        delta_e     = spa.atmospheric_refraction_correction(
            pressure[i], temp[i], e0, atmos_refract,
        )
        zenith[i] = 90.0-spa.topocentric_elevation_angle(e0, delta_e)

    return solardist, zenith
