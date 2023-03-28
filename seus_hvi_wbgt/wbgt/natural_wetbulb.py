"""
Various algorithms for computing natural wetbulb

"""

from .calc import relative_humidity

def malchaire(temp_a, temp_d, temp_w, temp_g ):
    """
    Compute natural wet bulb temperature

    An algorithm developed by Malchaire (1976) is used to compute the
    natural wet bulb temperature using the wet bulb, globe, ambient, and
    dew point temperature. Relative humidty is computerd using the
    ambient and dew point temperatures.

    Arguments:
        temp_a (ndarray) : Ambient temp (C)
        temp_d (ndarray) : Dew point temperature (C)
        temp_w (ndarray) : Wet bulb (C)
        temp_g (ndarray) : Globe temp (C)

    Reference:
        Malchaire, J. B., EVALUATION OF NATURAL WET BULB AND WET GLOBE THERMOMETERS.
            The Annals of Occupational Hygiene, Volume 19, Issue 3-4, December 1976,
            Pages 251–258, https://doi.org/10.1093/annhyg/19.3-4.251

    """

    relhum = relative_humidity( temp_a, temp_d )
    return (
        (0.16*(temp_g-temp_a) + 0.8)/200.0 *
        (560.0 - 2.0*relhum - 5.0*temp_a) - 0.8 + temp_w
    )

def hunter_minyard(temp_w, solar, speed):
    """
    Compute natural wet bulb temperature

    An algorithm shown in Hunter and Minyard (1999) to compute the 
    natural wetbulb temperture using the psychrometric wetbulb, 
    solar irradiance, and wind speed.
 
    Arguments:
        temp_w (ndarray) : Wet bulb globe temperature; degree C
        solar (ndarray) : Solar irradiance; W/m**2
        speed (ndarray) : Wind speed; m/s

    Reference:
        Hunter, Charles H., and C. Olivia Minyard. 
            Estimating wet bulb globe temperature using standard meteorological
            measurements." Proceedings of the Conference: 2nd Conference on
            Environmental Applications, Long Beach, CA, USA. Vol. 18. 1999.
 
    """

    #return temp_w + 0.021*solar - 0.42*speed + 1.93
    # Adjustment made based on the formula in the Boyer paper; see nws_boyer()
    return temp_w + 0.0021*solar - 0.43*speed + 1.93
def nws_boyer( temp_a, temp_w, solar, speed ):
    """
    Compute natural wet bulb temperature

    An algorithm shown in National Weather Service (NWS) documention
    from Broyer to compute the natural wetbulb temperture using
    psychrometric wetbulb, webulb depression,
    solar irradiance, and wind speed.
 
    Arguments:
        temp_a (ndarray) : Wet bulb globe temperature; degree C
        temp_w (ndarray) : Wet bulb globe temperature; degree C
        solar (ndarray) : Solar irradiance; W/m**2
        speed (ndarray) : Wind speed; m/s

    Reference:
        Boyer, Timothy R.
            NDFD Wet Bulb Globe Temperature Algorithm and Software Design
            https://vlab.noaa.gov/documents/6609493/7858379/
            NDFD+WBGT+Description+Document.pdf/fb89cc3a-0536-111f-f124-e4c93c746ef7?t=1642792547129 

    """

    return temp_w + 0.001651*solar - 0.09555*speed + 0.13235*(temp_a-temp_w) + 0.20249
