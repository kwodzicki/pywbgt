"""
Algorithms for computing psychrometric wetbulb temperature

"""

from cython.parallel import prange
import numpy

from libc.math cimport exp, fabsf

cimport numpy
cimport cython

from metpy.calc import relative_humidity_from_dewpoint as relative_humidity

cdef:
    float C0  =   26.66082
    float C1  =    0.0091379024
    float C2  = 6106.396

    float Cp  = 1004.0    # J/K/kg
    float L   =    2.54e6 # J/kg
    float eps =    0.622
    float f   = Cp/L/eps

cdef float saturation( float T ) nogil:

    return <float>exp( <float>C0 - C1*T - C2/T )

cdef float _iribarne_wb( float temp_a, float temp_d, float pres, int maxfev, float xtol) nogil:
    """
    Iteratively compute the wet-bulb temperature 

    Uses the Iribarne and Godson (1981) method to calculate the
    wet-bulb temperature given dry-bulb temperature, dew point
    temperature, and atmospheric pressure.

    Arguments:
        temp_a : Dry-bulb temperature (Kelvin)
        temp_d : Dew-point temperature (Kelvin)
        pres   : Atmoshperic pressure (hPa)
        maxfev : Maximum number of iterations to perform
        xtol : Calculation will terminate if the relative error between
            two consecutive iterates is at most xtol

    Returns:
        float : Wet-bulb temperature in Kelvin

    """

    cdef float e_sat, e_atm, e_wb, e_adj, temp_w, denom, numer, adjust 

    e_sat = saturation( temp_a )
    e_atm = saturation( temp_d )
    e_adj = (e_sat-e_atm)/(temp_a-temp_d)

    temp_w  = (temp_a*f*pres + temp_d*e_adj)/(f*pres + e_adj)

    while maxfev > 0:                                                             # While have NOT reached max iteration
        e_wb    = saturation( temp_w )                                                     # Compute saturation vapor pressure for wet-bulb temperature
        numer   = f*pres*(temp_a-temp_w) - (e_wb-e_atm)                                                # Computer error
        denom   = e_wb*(C1-C2/temp_w**2) - f*pres
        adjust  = numer / denom
        temp_w -= adjust 
        if fabsf(adjust) < xtol: return temp_w
        maxfev -= 1

    return 0.0/0.0

@cython.binding(True)
def stull( temp_a, temp_d ):
    """
    Wet bulb temperature from Stull method

    This formula for wet bulb temperature appears in:
        Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and 
            Air Temperature, Journal of Applied Meteorology and Climatology, 
            50(11), 2267-2269. Retrieved Jul 20, 2022, from 
            https://journals.ametsoc.org/view/journals/apme/50/11/jamc-d-11-0143.1.xml

    Arguments:
        temp_a (pint.Quantity) : Ambient (dry bulb) temperature 
        temp_d (pint.Quantity) : Dew point temperature

    """

    relhum = relative_humidity( temp_a, temp_d ).to('percent').magnitude
    temp_a = temp_a.to('degC').magnitude

    return (
        temp_a*numpy.arctan( 0.151977*(relhum + 8.313659)**(1.0/2.0) ) +
        numpy.arctan( temp_a + relhum ) - numpy.arctan( relhum - 1.676331 ) +
        0.00391838*relhum**(3.0/2.0)*numpy.arctan( 0.023101*relhum ) -
        4.686035
    )

@cython.binding(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def iribarne(temp_a, temp_d, pres=None, **kwargs):
    """
    Wet bulb temperature from Iribarne and Godson method

    Uses the Iribarne and Godson (1981) method to 
    Iteratively calculate the wet-bulb temperature given dry-bulb 
    temperature, dew point temperature, and atmospheric pressure.
 
    Arguments:
        temp_a (pint.Quantity): Dry-bulb temperature 
        temp_d (pint.Quantity) : Dew-point temperature

    Keyword arguments:
        pres  : Atmoshperic pressure (hPa)
        maxfev : Maximum number of iterations to perform
        xtol : Calculation will terminate if the relative error between
            two consecutive iterates is at most xtol

    Returns:
        numpy.ndarray : Wet bulb temperature(s)

    """

    cdef:
        Py_ssize_t i, size = temp_a.size
        int maxfev = kwargs.get('maxfev', 25)
        float xtol = kwargs.get('xtol',   0.02)

    out = numpy.empty( size, dtype = numpy.float32 )
    cdef:
        float [::1] out_view    = out
        float [::1] temp_a_view = (
            temp_a
            .to('kelvin')
            .magnitude
            .astype( numpy.float32 )
        )
        float [::1] temp_d_view = (
            temp_d
            .to('kelvin')
            .magnitude
            .astype( numpy.float32 )
        )
        float [::1] pres_view   = (
            numpy.full( size, 1013.25, dtype=numpy.float32 )
            if pres is None else
            pres.astype(  numpy.float32 )
        )

    for i in prange( size, nogil=True ):
        out_view[i] = _iribarne_wb(
            temp_a_view[i],
            temp_d_view[i],
            pres_view[i],
            maxfev,
            xtol,
        ) - 273.15

    return out
