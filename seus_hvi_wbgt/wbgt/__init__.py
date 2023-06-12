"""
Package for estimating wetbulb globe temperature

Various algorithms for estimating wetbulb globe
temperature from standard meteorological variables.
 
"""

from .liljegren import wetbulb_globe as liljegrenWBGT
from .bernard   import wetbulb_globe as bernardWBGT
from .dimiceli  import wetbulb_globe as dimiceliWBGT

def wbgt( method, *args, **kwargs ):
    """
    Estimate wet bulb globe temperature

    Wrapper for the various WBGT algorithms provided
    by this package. Set the method to use for 
    estimating WBGT and you're off

    Arguments:
        method (str) : name of the method to use
        lat (ndarray) : Latitude corresponding to data values (decimal).
            Can be one (1) element array; will be expanded to match dates/data
        lon (ndarray) : Longitude correspondning to data values (decimal).
            Can be one (1) element array; will be expanded to match dates/data
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        solar (Quantity) : solar irradiance; units of any power over area
        pres (Qantity) : barometric pressure; units of pressure
        temp_air (Quantity) : air (dry bulb) temperature; units of temperature
        temp_dew (Quantity) : Dew point temperature; units of temperature
        speed (Quatity) : wind speed; units of speed

    Keyword arguments:
        f_db (float) : Direct beam radiation from the sun. Valid for
            the Dimiceli and Bernard methods. Type: fraction
        cosz (float) : Cosine of solar zenith angle. Valid for 
            Dimiceli and Bernard methods
        zspeed (Quantity) : Height of the wind speed measurment.
            Default is 10 meters
        wetbulb (str) : Name of wet bulb algorithm to use in the Dimiceli
            algorithm. Valid options are:
            {dimiceli, stull} DEFAULT = dimiceli
        natural_wetbulb (str) : Name of the natural wet bulb algorithm to use
            in the Dimiceli algorithm. Valid options are:
            {malchaire, hunter_minyard} DEFAULT = malchaire
        urban (ndarray) : Boolean flag indicating if "urban" (1) or
            "rural" (0) for wind speed power law exponent. Valid for the
            Liljegren algorithm.
            Can be one (1) element array; will be expanded to match dates/data
        gmt (ndarray) LST-GMT difference  (hours; negative in USA). Valid for
            the Liljegren algorithm.
        avg (ndarray) : averaging time of the meteorological inputs (minutes).
            Valid for the Liljegren algorithm.
        dT (Quantity) : Vertical temperature difference; upper minus lower;
            unit of temperature. Valid for the Liljegren algorithm.
        use_spa (bool) : If set, use the National Renewable Energy
            Laboratory (NREL) Solar Position Algorithm (SPA) to determine
            sun position. Default is to use the build-it, low precision model.
        d_globe (Quantity) : Diameter of the black globe thermometer
            unit of distance. Valid for the Liljegren algorithm.

    Returns:
        dict :
            - Tg : Globe temperatures as Quantity
            - Tpsy : psychrometric wet bulb temperatures as Quantity
            - Tnwb : Natural wet bulb temperatures as Quantity
            - Twbg : Wet bulb-globe temperatures as Quantity
            - solar : Adjusted solar irradiance as Quantity.
            - speed : Estimated 2m wind speed as Quantity:
                will be same as input if already 2m t

    """

    method = method.lower()
    if method == 'liljegren':
        return liljegrenWBGT( *args, **kwargs )
    if method == 'bernard':
        return bernardWBGT( *args, **kwargs )
    if method == 'dimiceli':
        return dimiceliWBGT( *args, **kwargs )

    raise Exception( f'Unsupported WBGT method : {method}' )
