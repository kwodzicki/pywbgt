"""
Solar position and direct beam fraction

Python implementation of the Liljegren solar position and direct beam
fraction calculations utilizing the Python pvlib implementation of the
National Renewable Energy Laboratory (NREL) Solar Position Algorithm (SPA).

"""

import numpy as np
from pvlib import spa

from .utils import datetime_adjust
from .liljegren import (
    LILJEGREN_CZA_MIN, LILJEGREN_NORMSOLAR_MAX, LILJEGREN_SOLAR_CONST,
)

# Default values for some parameters
ELEV = 0.00
PRESSURE = 1013.25
TEMP = 15.00


def solar_parameters(
    datetime,
    lat,
    lon,
    solar,
    gmt=None,
    avg=None,
    elev=None,
    pressure=None,
    temp=None,
    **kwargs,
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
        lat (ndarray) : Latitude of location(s) to compute parameters for;
            decimal
        lon (ndarray) : Longitude of location(s) to compute parameters for;
            decimal
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
    ) / 1.0e9

    ntime = datetime.shape[0]

    # If input latitude is only one (1) element, assume lon and urban are
    # also one (1) element and expand all to match size of data
    if lat.size <= 1:
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

    return _solar_parameters(
        datetime,
        lat,
        lon,
        solar,
        elev,
        pressure,
        temp,
        0,
        0,
    )


def _solar_parameters(
    unixtime,
    lat,
    lon,
    solar,
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

    jd = spa.julian_day(unixtime)
    jde = spa.julian_ephemeris_day(jd, delta_t)
    jc = spa.julian_century(jd)
    jce = spa.julian_ephemeris_century(jde)
    jme = spa.julian_ephemeris_millennium(jce)
    R = heliocentric_radius_vector(jme)

    L = heliocentric_longitude(jme)
    B = heliocentric_latitude(jme)
    Theta = spa.geocentric_longitude(L)
    beta = spa.geocentric_latitude(B)
    delta_psi, delta_epsilon = longitude_obliquity_nutation(
        jce,
        spa.mean_elongation(jce),
        spa.mean_anomaly_sun(jce),
        spa.mean_anomaly_moon(jce),
        spa.moon_argument_latitude(jce),
        spa.moon_ascending_longitude(jce),
    )

    epsilon = spa.true_ecliptic_obliquity(
        spa.mean_ecliptic_obliquity(jme),
        delta_epsilon,
    )
    delta_tau = spa.aberration_correction(R)

    lamd = spa.apparent_sun_longitude(Theta, delta_psi, delta_tau)
    v = spa.apparent_sidereal_time(
        spa.mean_sidereal_time(jd, jc),
        delta_psi,
        epsilon,
    )

    alpha = spa.geocentric_sun_right_ascension(lamd, epsilon, beta)
    delta = spa.geocentric_sun_declination(lamd, epsilon, beta)

    H = spa.local_hour_angle(v, lon, alpha)
    xi = spa.equatorial_horizontal_parallax(R)
    u = spa.uterm(lat)
    x = spa.xterm(u, lat, elev)
    y = spa.yterm(u, lat, elev)

    delta_alpha = spa.parallax_sun_right_ascension(x, xi, H, delta)
    delta_prime = spa.topocentric_sun_declination(
        delta,
        x,
        y,
        xi,
        delta_alpha,
        H,
    )
    e0 = spa.topocentric_elevation_angle_without_atmosphere(
        lat,
        delta_prime,
        spa.topocentric_local_hour_angle(H, delta_alpha),
    )
    delta_e = spa.atmospheric_refraction_correction(
        pressure, temp, e0, atmos_refract,
    )

    cza = np.cos(
        np.deg2rad(90.0 - spa.topocentric_elevation_angle(e0, delta_e))
    )

    fdir = np.zeros(cza.shape)
    idx = cza >= LILJEGREN_CZA_MIN
    if idx.any():
        toasolar = LILJEGREN_SOLAR_CONST * cza[idx].clip(min=0.0) / R[idx]**2

        # Limit maximum value of norm solar
        normsolar = (solar[idx] / toasolar).clip(max=LILJEGREN_NORMSOLAR_MAX)

        solar[idx] = normsolar * toasolar
        iidx = normsolar > 0.0
        if iidx.any():
            _fdir = np.zeros(normsolar.shape)
            _fdir[iidx] = np.exp(
                3.0 - 1.34 * normsolar[iidx] - 1.65 / normsolar[iidx]
            )
            fdir[idx] = _fdir

    solar[~idx] = 0
    return solar, cza, fdir.clip(min=0.0, max=0.9)


def sum_mult_cos_add_mult(arr, x):
    """From pvlib.spa, updates for array operations"""

    nn = (1,) * x.ndim + arr.shape[:1]
    return (
        arr[:, 0].reshape(nn)
        * np.cos(
            arr[:, 1].reshape(nn)
            + arr[:, 2].reshape(nn) * x.reshape(x.shape + (1,))
        )
    ).sum(axis=-1)


def heliocentric_radius_vector(jme):
    """From pvlib.spa, updates for array operations"""

    r0 = sum_mult_cos_add_mult(spa.R0, jme)
    r1 = sum_mult_cos_add_mult(spa.R1, jme)
    r2 = sum_mult_cos_add_mult(spa.R2, jme)
    r3 = sum_mult_cos_add_mult(spa.R3, jme)
    r4 = sum_mult_cos_add_mult(spa.R4, jme)

    return (r0 + r1 * jme + r2 * jme**2 + r3 * jme**3 + r4 * jme**4) / 10**8


def heliocentric_longitude(jme):
    """From pvlib.spa, updates for array operations"""

    l0 = sum_mult_cos_add_mult(spa.L0, jme)
    l1 = sum_mult_cos_add_mult(spa.L1, jme)
    l2 = sum_mult_cos_add_mult(spa.L2, jme)
    l3 = sum_mult_cos_add_mult(spa.L3, jme)
    l4 = sum_mult_cos_add_mult(spa.L4, jme)
    l5 = sum_mult_cos_add_mult(spa.L5, jme)

    l_rad = (
        l0 + l1 * jme + l2 * jme**2 + l3 * jme**3 + l4 * jme**4
        + l5 * jme**5
    ) / 10**8
    return np.rad2deg(l_rad) % 360


def heliocentric_latitude(jme):
    """From pvlib.spa, updates for array operations"""

    b0 = sum_mult_cos_add_mult(spa.B0, jme)
    b1 = sum_mult_cos_add_mult(spa.B1, jme)

    b_rad = (b0 + b1 * jme) / 10**8
    return np.rad2deg(b_rad)


def longitude_obliquity_nutation(
    julian_ephemeris_century,
    x0,
    x1,
    x2,
    x3,
    x4,
):
    """From pvlib.spa, updates for array operations"""

    nn = (
        (1,) * julian_ephemeris_century.ndim
        + spa.NUTATION_YTERM_ARRAY.shape[:1]
    )
    mm = julian_ephemeris_century.shape + (1,)

    julian_ephemeris_century = julian_ephemeris_century.reshape(mm)

    a = spa.NUTATION_ABCD_ARRAY[:, 0].reshape(nn)
    b = spa.NUTATION_ABCD_ARRAY[:, 1].reshape(nn)
    c = spa.NUTATION_ABCD_ARRAY[:, 2].reshape(nn)
    d = spa.NUTATION_ABCD_ARRAY[:, 3].reshape(nn)

    arg = np.radians(
        spa.NUTATION_YTERM_ARRAY[:, 0].reshape(nn) * x0.reshape(mm)
        + spa.NUTATION_YTERM_ARRAY[:, 1].reshape(nn) * x1.reshape(mm)
        + spa.NUTATION_YTERM_ARRAY[:, 2].reshape(nn) * x2.reshape(mm)
        + spa.NUTATION_YTERM_ARRAY[:, 3].reshape(nn) * x3.reshape(mm)
        + spa.NUTATION_YTERM_ARRAY[:, 4].reshape(nn) * x4.reshape(mm)
    )
    delta_psi_sum = (
        (a + b * julian_ephemeris_century) * np.sin(arg)
    ).sum(axis=-1)

    delta_eps_sum = (
        (c + d * julian_ephemeris_century) * np.cos(arg)
    ).sum(axis=-1)

    return delta_psi_sum / 36000000, delta_eps_sum / 36000000
