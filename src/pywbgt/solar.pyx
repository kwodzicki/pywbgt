"""
Solar position and direct beam fraction

Python implementation of the Liljegren solar position and direct beam
fraction calculations utilizing the Python pvlib implementation of the
National Renewable Energy Laboratory (NREL) Solar Position Algorithm (SPA).

"""

import numpy as np
cimport numpy as cnp

from libc.math cimport cos, sin, fmod, M_PI
import cython
from cython.parallel cimport prange

cnp.import_array()

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

    assert datetime.ndim < 2
    assert lat.ndim < 2
    assert lon.ndim < 2
    assert solar.ndim < 2

    datetime = (
        datetime_adjust(datetime, gmt, avg)
        .astype(np.int64)
    ) / 1.0e9

    ntime = datetime.shape

    # If input latitude is only one (1) element, assume lon and urban are
    # also one (1) element and expand all to match size of data
    if lat.size <= 1:
        lat = np.full(ntime, lat)

    if lon.size <= 1:
        lon = np.full(ntime, lon)

    if elev is None:
        elev = np.full(ntime, ELEV)
    elif elev.size <= 1:
        elev = np.full(ntime, elev)

    if pressure is None:
        pressure = np.full(ntime, PRESSURE)
    elif pressure.size <= 1:
        pressure = np.full(ntime, pressure)

    if temp is None:
        temp = np.full(ntime, TEMP)
    elif temp.size <= 1:
        temp = np.full(ntime, temp)

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
    jce = spa.julian_ephemeris_century(
        spa.julian_ephemeris_day(jd, delta_t)
    )
    jme = spa.julian_ephemeris_millennium(jce)

    R = heliocentric_radius_vector(jme)

    Theta = spa.geocentric_longitude(
        heliocentric_longitude(jme)
    )
    beta = spa.geocentric_latitude(
        heliocentric_latitude(jme)
    )
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
        spa.mean_sidereal_time(jd, spa.julian_century(jd)),
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
            _fdir = np.zeros_like(normsolar)
            _fdir[iidx] = np.exp(
                3.0 - 1.34 * normsolar[iidx] - 1.65 / normsolar[iidx]
            )
            fdir[idx] = _fdir

    solar[~idx] = 0
    return solar, cza, fdir.clip(min=0.0, max=0.9)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def heliocentric_radius_vector(double [::1] jme):
    """From pvlib.spa, updates for array operations"""


    cdef Py_ssize_t i    
    cdef Py_ssize_t n_jme = jme.size
    cdef Py_ssize_t n_R0 = spa.R0.shape[0]
    cdef Py_ssize_t n_R1 = spa.R1.shape[0]
    cdef Py_ssize_t n_R2 = spa.R2.shape[0]
    cdef Py_ssize_t n_R3 = spa.R3.shape[0]
    cdef Py_ssize_t n_R4 = spa.R4.shape[0]

    cdef double jme1, jme2, jme3, jme4

    cdef double [:, ::1] R0 = spa.R0
    cdef double [:, ::1] R1 = spa.R1
    cdef double [:, ::1] R2 = spa.R2
    cdef double [:, ::1] R3 = spa.R3
    cdef double [:, ::1] R4 = spa.R4
    
    cdef double [::1] res = np.zeros(n_jme, dtype=np.float64)

    with nogil:
        for i in prange(n_jme):
            jme1 = jme[i]
            jme2 = jme1 * jme1
            jme3 = jme2 * jme1
            jme4 = jme3 * jme1

            res[i] = csum_mult_cos_add_mult(n_R0, R0, jme1)
            res[i] += csum_mult_cos_add_mult(n_R1, R1, jme1) * jme1
            res[i] += csum_mult_cos_add_mult(n_R2, R2, jme1) * jme2
            res[i] += csum_mult_cos_add_mult(n_R3, R3, jme1) * jme3
            res[i] += csum_mult_cos_add_mult(n_R4, R4, jme1) * jme4
            res[i] *= 1.0e-8
 
    return np.asarray(res)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def heliocentric_longitude(double [::1] jme):
    """From pvlib.spa, updates for array operations"""


    cdef Py_ssize_t i
    cdef Py_ssize_t n_jme = jme.size
    cdef Py_ssize_t n_L0 = spa.L0.shape[0]
    cdef Py_ssize_t n_L1 = spa.L1.shape[0]
    cdef Py_ssize_t n_L2 = spa.L2.shape[0]
    cdef Py_ssize_t n_L3 = spa.L3.shape[0]
    cdef Py_ssize_t n_L4 = spa.L4.shape[0]
    cdef Py_ssize_t n_L5 = spa.L5.shape[0]

    cdef double jme1, jme2, jme3, jme4, jme5

    cdef double [:, ::1] L0 = spa.L0
    cdef double [:, ::1] L1 = spa.L1
    cdef double [:, ::1] L2 = spa.L2
    cdef double [:, ::1] L3 = spa.L3
    cdef double [:, ::1] L4 = spa.L4
    cdef double [:, ::1] L5 = spa.L5
    
    cdef double [::1] res = np.zeros(n_jme, dtype=np.float64)

    with nogil:
        for i in prange(n_jme):
            jme1 = jme[i]
            jme2 = jme1 * jme1
            jme3 = jme2 * jme1
            jme4 = jme3 * jme1
            jme5 = jme4 * jme1

            res[i] = csum_mult_cos_add_mult(n_L0, L0, jme1)
            res[i] += csum_mult_cos_add_mult(n_L1, L1, jme1) * jme1
            res[i] += csum_mult_cos_add_mult(n_L2, L2, jme1) * jme2
            res[i] += csum_mult_cos_add_mult(n_L3, L3, jme1) * jme3
            res[i] += csum_mult_cos_add_mult(n_L4, L4, jme1) * jme4
            res[i] += csum_mult_cos_add_mult(n_L5, L5, jme1) * jme5
            res[i] = pymod((res[i] * 1.0e-8) * 180 / M_PI, 360.0)
 
    return np.asarray(res)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def heliocentric_latitude(double [::1] jme):
    """From pvlib.spa, updates for array operations"""

    cdef Py_ssize_t i
    cdef Py_ssize_t n_jme = jme.size
    cdef Py_ssize_t n_B0 = spa.B0.shape[0]
    cdef Py_ssize_t n_B1 = spa.B1.shape[0]
    
    cdef double [::1] res = np.zeros(n_jme, dtype=np.float64)

    cdef double [:, ::1] B0 = spa.B0
    cdef double [:, ::1] B1 = spa.B1

    with nogil:
        for i in prange(n_jme):
            res[i] = csum_mult_cos_add_mult(n_B0, B0, jme[i])
            res[i] += csum_mult_cos_add_mult(n_B1, B1, jme[i]) * jme[i]
            res[i] = (res[i] * 1.0e-8) * 180 / M_PI

    return np.asarray(res)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def longitude_obliquity_nutation(
    double[::1] jec,
    double[::1] x0,
    double[::1] x1,
    double[::1] x2,
    double[::1] x3,
    double[::1] x4,
):
    """From pvlib.spa, updates for array operations"""


    cdef Py_ssize_t i, j
    cdef Py_ssize_t n_arr = spa.NUTATION_YTERM_ARRAY.shape[0]
    cdef Py_ssize_t n_x = jec.shape[0]
    cdef double factor = 1.0 / 36000000
    cdef double radians = M_PI / 180.0
    cdef double arg
    cdef double[:, ::1] abcd = spa.NUTATION_ABCD_ARRAY
    cdef long[:, ::1] yterm = spa.NUTATION_YTERM_ARRAY

    cdef cnp.ndarray[cnp.float64_t, ndim=1] delta_psi = np.zeros(
        n_x, dtype=np.float64,
    )
    cdef cnp.ndarray[cnp.float64_t, ndim=1] delta_eps = np.zeros(
        n_x, dtype=np.float64,
    )

    # Parallel loop over each element in x
    for i in prange(n_x, nogil=True):
        arg = 0.0
        for j in range(n_arr):
            arg = radians * ( 
                yterm[j, 0] * x0[i]
                + yterm[j, 1] * x1[i]
                + yterm[j, 2] * x2[i]
                + yterm[j, 3] * x3[i]
                + yterm[j, 4] * x4[i]
            )
            delta_psi[i] += (abcd[j, 0] + abcd[j, 1] * jec[i]) * sin(arg)
            delta_eps[i] += (abcd[j, 2] + abcd[j, 3] * jec[i]) * cos(arg)

        delta_psi[i] *= factor
        delta_eps[i] *= factor

    return delta_psi, delta_eps


@cython.boundscheck(False)
cdef double csum_mult_cos_add_mult(
    Py_ssize_t n_arr,
    double [:, ::1] arr,
    double x
) nogil:
    cdef Py_ssize_t i
    cdef double res = 0.0

    # Parallel loop over each element in x
    for i in range(n_arr):
        res += arr[i, 0] * cos(arr[i, 1] + arr[i, 2] * x)

    return res


@cython.cdivision(True)
cdef inline double pymod(double a, double b) nogil:
    """Python-style floating point modulus"""
    cdef double r = fmod(a, b)
    if (r != 0.0) and ((r < 0.0) != (b < 0.0)):
        r += b
    return r
