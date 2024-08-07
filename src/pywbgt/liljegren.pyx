from cython.parallel import prange
import numpy

#from libc.stdio cimport printf
from libc.math cimport fmaxf

cimport numpy
cimport cython

from metpy.units import units
from metpy.calc import relative_humidity_from_dewpoint as rhTd

from .cliljegren cimport *

# Expose cython constants to python
LILJEGREN_D_GLOBE       = units.Quantity(_D_GLOBE,   'meter')
LILJEGREN_MIN_SPEED     = units.Quantity(_MIN_SPEED, 'meter/second')
LILJEGREN_CZA_MIN       = _CZA_MIN
LILJEGREN_NORMSOLAR_MAX = _NORMSOLAR_MAX
LILJEGREN_SOLAR_CONST   = _SOLAR_CONST

from .constants import MIN_SPEED
from .solar import solar_parameters as sparms

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.always_allow_keywords(True)
def conv_heat_trans_coeff( temp_air, pres, speed, float diameter=_D_GLOBE ):
    """
    Compute convective heat transfer coefficient 
  
    Wrapper for the h_sphere_in_air C function from WBGT v1.1
  
    Arguments:
        temp_air (ndarray) : Ambient air temperature in Kelvin
        pres (ndarray) : Barometric pressure in hPa/mb
        speed (ndarray) : Air/wind seed in m/s
  
    Keyword arguments:
        diameter (float) : Diameter of the sphere in meters; default value
            is 0.0508 and is taken from the C code
  
    Returns:
        ndarray : convective heat transfer coefficients

    """
  
    cdef Py_ssize_t i, size = temp_air.shape[0]
  
    h = numpy.empty( size, dtype = numpy.float32 )
  
    cdef:
        float [::1] hView        = h # Initialize array to write data to
        float [::1] temp_airView = temp_air.astype( numpy.float32 ) 
        float [::1] presView     = pres.astype( numpy.float32 )
        float [::1] speedView    = speed.astype( numpy.float32 )
  
    for i in prange( size, nogil=True ):                                          # Iterate over all values in parallel
        hView[i] = h_sphere_in_air( 
            diameter, temp_airView[i], presView[i], speedView[i]
        )
  
    return h                                                     # Reshape to same shape as temp_air 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
@cython.always_allow_keywords(True)
def globe_temperature(
        temp_air, temp_dew, pres, speed, solar, fdir, cza,
        d_globe = None,
    ):
    """
    Compute globe temperature using Liljegren method

    Arguments:
        temp_air (ndarray) : Air (dry bulb) temperature; degree Celsius
        temp_dew (ndarray) : Dew point temperature; degree Celsius
        pres (ndarray) : Barometric pressure; hPa
        speed (ndarray) : wind speed, m/s
        solar (ndarray) : Solar irradiance, W/m**2
        fdir (ndarray) : Fraction of solar irradiance due to direct beam
        cza (ndarray) : Cosine of solar zenith angle

    Returns:
        ndarray : Globe temperature

    """

    cdef:
        float _d_globe
        float [::1] temp_airView = (temp_air + 273.15).astype(  numpy.float32 )
        float [::1] presView     = pres.astype(  numpy.float32 )
        float [::1] speedView    = speed.astype( numpy.float32 )
        float [::1] solarView    = solar.astype( numpy.float32 )
        float [::1] fdirView     = fdir.astype(  numpy.float32 )
        float [::1] czaView      = cza.astype(   numpy.float32 )
        float [::1] relhumView   = (
            rhTd(
                units.Quantity(temp_air, 'degC'),
                units.Quantity(temp_dew, 'degC'),
            )
            .magnitude
            .astype( numpy.float32 )
        )

        Py_ssize_t i, size = temp_air.shape[0]

    if d_globe is None:
        _d_globe = _D_GLOBE                                                          # Use default value from C source code
    else:
        _d_globe = d_globe.to('meter').astype(numpy.float32).magnitude              # Ensure is in units of meter, convert to 32-bit float, and get magnitude

     
    out = numpy.empty( size, dtype=numpy.float32 )
    cdef float [:] outView = out

    for i in prange( size, nogil=True ):
        outView[i] = Tglobe(
            temp_airView[i], relhumView[i], presView[i], 
            speedView[i], solarView[i], fdirView[i],
            czaView[i], _d_globe,
        )

    return out 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
@cython.always_allow_keywords(True)
def psychrometric_wetbulb( temp_air, temp_dew, pres ):
    """
    Compute psychrometeric wet bulb temperature using Liljegren method

    Arguments:
        temp_air (ndarray) : Air (dry bulb) temperature; degree Celsius
        temp_dew (Quantity) : Dew point temperature; units of temperature
        pres (ndarray) : Barometric pressure; hPa

    Returns:
        ndarray : pyschrometric Wet bulb temperature

    """

    cdef:
        float [::1] temp_airView = (temp_air + 273.15).astype(  numpy.float32 )
        float [::1] presView     = pres.astype(  numpy.float32 )
        float [::1] relhumView   = (
            rhTd(
                units.Quantity(temp_air, 'degC'),
                units.Quantity(temp_dew, 'degC'),
            )
            .magnitude
            .astype( numpy.float32 )
        )

        float tmp, fill = 0.0
        int   rad  = 0

        Py_ssize_t i, size = temp_air.shape[0]

    out = numpy.full( size, numpy.nan, dtype=numpy.float32 )
    cdef float [::1] outView = out

    for i in prange( size, nogil=True ):
        tmp = Twb(
            temp_airView[i], 
            relhumView[i], 
            presView[i], 
            fill, fill, fill, fill, rad
        )
        if tmp > -9999.0:
            outView[i] = tmp

    return out 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
@cython.always_allow_keywords(True)
def natural_wetbulb( temp_air, temp_dew, pres, speed, solar, fdir, cza ):
    """
    Compute natural wet bulb temperature using Liljegren method

    Arguments:
      temp_air (ndarray) : Air (dry bulb) temperature; degree Celsius
      temp_dew (ndarray) : Dew point temperature; units of temperature
      pres (ndarray) : Barometric pressure; hPa
      speed (ndarray) : wind speed, m/s
      solar (ndarray) : Solar irradiance, W/m**2
      fdir (ndarray) : Fraction of solar irradiance due to direct beam
      cza (ndarray) : Cosine of solar zenith angle

    Returns:
      ndarray : natural Wet bulb temperature

    """

    cdef:
        float [::1] temp_airView = (temp_air + 273.15).astype(  numpy.float32 )
        float [::1] presView     = pres.astype(  numpy.float32 )
        float [::1] speedView    = speed.astype( numpy.float32 )
        float [::1] solarView    = solar.astype( numpy.float32 )
        float [::1] fdirView     = fdir.astype(  numpy.float32 )
        float [::1] czaView      = cza.astype(   numpy.float32 )
        float [::1] relhumView   = (
            rhTd(
                units.Quantity(temp_air, 'degC'),
                units.Quantity(temp_dew, 'degC'),
            )
            .magnitude
            .astype( numpy.float32 )
        )

        float tmp
        int rad = 1

        Py_ssize_t i, size = temp_air.shape[0]

    out = numpy.full( size, numpy.nan, dtype=numpy.float32 )
    cdef float [:] outView = out

    for i in prange( size, nogil=True ):
        tmp = Twb( 
            temp_airView[i], relhumView[i], presView[i], 
            speedView[i], solarView[i], fdirView[i], czaView[i], rad
        )
        if tmp > -9999.0:
            outView[i] = tmp

    return out 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
@cython.always_allow_keywords(True)
def wetbulb_globe(
        datetime, lat, lon,
        solar, pres, temp_air, temp_dew, speed,
        urban     = None, 
        gmt       = None,
        avg       = None,
        zspeed    = None,
        dT        = None,
        min_speed = None,
        d_globe   = None,
        **kwargs,
    ): 

    """
    Calculate the outdoor wet bulb-globe temperature

    Cython wrapper for Liljegren C code for calculating the outdoor
    wet bulb-globe temperature,, which is the weighted sum of the\
    air temperature (dry bulb), the globe temperature, and the 
    natural wet bulb temperature:

        Twbg = 0.1 * temp_air + 0.7 * Tnwb + 0.2 * Tg.

    The program predicts Tnwb and Tg using meteorological input data
    then combines the results to produce Twbg.

    Arguments:
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        lat (ndarray) : Latitude corresponding to data values (decimal).
            Can be one (1) element array; will be expanded to match dates/data
        lon (ndarray) : Longitude correspondning to data values (decimal).
            Can be one (1) element array; will be expanded to match dates/data
        solar (Quantity) : solar irradiance; units of any power over area
        pres (Qantity) : barometric pressure; units of pressure
        temp_air (Quantity) : air (dry bulb) temperature; units of temperature
        temp_dew (Quantity) : Dew point temperature; units of temperature
        speed (Quatity) : wind speed; units of speed

    Keyword arguments:
        urban (ndarray) : Boolean flag indicating if "urban" (1) or
            "rural" (0) for wind speed power law exponent
            Can be one (1) element array; will be expanded to match dates/data
        gmt (ndarray) LST-GMT difference  (hours; negative in USA)
        avg (ndarray) : averaging time of the meteorological inputs (minutes)
        zspeed (Quantity) : height of wind speed measurement; unit of distance
        dT (Quantity) : Vertical temperature difference; upper minus lower;
            unit of temperature
        min_speed (Quantity) : Sets the minimum speed for the height-adjusted
            wind speed. If this keyword is set, the larger of input value and
            LILJEGREN_MIN_SPEED is used. The default value is MIN_SPEED, which
            is 2 knots, and is larger than LILJEGREN_MIN_SPEED.
        d_globe (Quantity) : Diameter of the black globe thermometer
            unit of distance

    Returns:
        dict :
            - Tg : Globe temperatures as Quantity
            - Tpsy : psychrometric wet bulb temperatures as Quantity
            - Tnwb : Natural wet bulb temperatures as Quantity
            - Twbg : Wet bulb-globe temperatures as Quantity
            - Speed : Estimated 2m wind speed as Quantity; will be same as input if already 2m wind speed 
            - min_speed : Minimum speed that adjusted wind speed is clipped to as Quantity 

    Reference: 
        Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp, and R. Sharp:
            Modeling the Wet Bulb Globe Temperature Using Standard Meteorological
            Measurements. The Journal of Occupational and Environmental Hygiene,
            vol. 5:10, pp. 645-655, 2008.

    """

    # Define size of output arrays based on size of input
    cdef:
        Py_ssize_t i, size = datetime.shape[0]
        float _min_speed, _d_globe

    # Generate or repeat urban value based on input
    if urban is None: 
        urban = numpy.zeros( size, dtype = numpy.int32 )
    elif len(urban) == 1:
        urban = urban.repeat( size )

    # Generature or repeat wind speed observation height based on input
    if zspeed is None:
        zspeed = (
            units.meter * 
            numpy.full(size, 10.0, dtype=numpy.float32) 
        )
    elif not hasattr( zspeed, 'size' ):
        zspeed = (
            zspeed.units *
            numpy.full(size, zspeed.magnitude, dtype=numpy.float32)
        )
    elif zspeed.size != size:
        raise Exception("Size mismatch between 'zspeed' and other variables!")

    # Set default temperature differential between 10m and 2m samples
    if dT is None:
        dT = (
            units.degree_Celsius * 
            numpy.full(size, -1.0, dtype=numpy.float32)
        ) 

    # Set default value for minimum speed
    if min_speed is None:
        # Greater of MIN_SPEED for package (2 knots) and absolute min speed for algorithm
        _min_speed = max(
            MIN_SPEED,
            LILJEGREN_MIN_SPEED
        ).to('meter/second').magnitude
    else:
        # Greater of user input min_speed and absolute min speed for algorithm
        _min_speed = max(
            min_speed,
            LILJEGREN_MIN_SPEED,
        ).to('meter/second').magnitude
    #printf("min speed %f\n", _min_speed)

    # Set default black globe thermometer diameter
    if d_globe is None:
        _d_globe = _D_GLOBE                                                          # Use default value from C source code
    else:
        _d_globe = d_globe.to('meter').astype(numpy.float32).magnitude              # Ensure is in units of meter, convert to 32-bit float, and get magnitude
    
    #printf("%f\n", _d_lobe)  
    if len( lat ) == 1:                                                         # If input latitude is only one (1) element, assume lon and urban are also one (1) element and expand all to match size of data
        lat = lat.repeat( size )
        lon = lon.repeat( size )

    #solar_adj, cza, fdir = solar_parameters(
    solar_adj, cza, fdir = sparms(
        datetime,
        lat,
        lon,
        solar.to('watt/m**2').magnitude,
        gmt,
        avg,
        **kwargs,
    )

    # Define output array for storing WBGT status code
    out = numpy.full( (6, size), numpy.nan, dtype = numpy.float32 )

    # Define array views for faster/parallel iteration
    cdef:
        int daytime, stability_class
        float Tg, Tpsy, Tnwb, Twbg, est_speed  
    
        float [:,::1] outView = out 

        float [::1] latView    = lat.astype(   numpy.float32 )
        float [::1] lonView    = lon.astype(   numpy.float32 )
        int   [::1] urbanView  = urban.astype( numpy.int32   )

        float [::1] solar_adjView   =  solar_adj.astype( numpy.float32 )
        float [::1] czaView         = cza.astype( numpy.float32)
        float [::1] fdirView        = fdir.astype( numpy.float32)

        float [::1] presView        =     pres.to('hPa'           ).magnitude.astype( numpy.float32 )
        float [::1] temp_airView    = temp_air.to('kelvin'        ).magnitude.astype( numpy.float32 )
        float [::1] speedView       =    speed.to('m/s'           ).magnitude.astype( numpy.float32 )
        float [::1] zspeedView      =   zspeed.to('meter'         ).magnitude.astype( numpy.float32 )
        float [::1] dTView          =       dT.to('degree_Celsius').magnitude.astype( numpy.float32 )
        float [::1] relhumView      = (
            rhTd(temp_air, temp_dew)
            .magnitude
            .astype( numpy.float32 )
        )


    # Iterate (in parallel) over all values in the input arrays
    for i in prange( size, nogil=True ):
        if zspeedView[i] == _REF_HEIGHT:
            est_speed = fmaxf(speedView[i], _min_speed)
        else:
            if czaView[i] > 0.0:
                daytime = 1
            stability_class = stab_srdt(
                daytime,
                speedView[i],
                solar_adjView[i],
                dTView[i],
            )
            est_speed = est_wind_speed(
                speedView[i],
                zspeedView[i],
                stability_class,
                urbanView[i],
                _min_speed,
            )
 
        # Ensure that variables are defined for each thread
        Tg = Tglobe(
            temp_airView[i],
            relhumView[i],
            presView[i],
            est_speed,
            solar_adjView[i],
            fdirView[i],
            czaView[i],
            _d_globe,
        )
        if Tg == -9999:
            continue

        Tnwb = Twb(
            temp_airView[i],
            relhumView[i],
            presView[i],
            est_speed,
            solar_adjView[i],
            fdirView[i],
            czaView[i],
            1,
        )
        if Tnwb == -9999:
            continue

        Tpsy = Twb(
            temp_airView[i],
            relhumView[i],
            presView[i],
            est_speed,
            solar_adjView[i],
            fdirView[i],
            czaView[i],
            0,
        )
 
        outView[0,i] = Tg
        outView[1,i] = Tpsy
        outView[2,i] = Tnwb
        outView[3,i] = 0.1*(temp_airView[i]-273.15) + 0.2*Tg + 0.7*Tnwb
        outView[4,i] = solar_adjView[i]
        outView[5,i] = est_speed

    # Return dict with unit-aware values
    return {
        'Tg'        : units.Quantity(out[0,:],   'degree_Celsius'),
        'Tpsy'      : units.Quantity(out[1,:],   'degree_Celsius'),
        'Tnwb'      : units.Quantity(out[2,:],   'degree_Celsius'),
        'Twbg'      : units.Quantity(out[3,:],   'degree_Celsius'),
        'solar'     : units.Quantity(out[4,:],   'watt/meter**2'),
        'speed'     : units.Quantity(out[5,:],   'meter/second'),
        'min_speed' : units.Quantity(_min_speed, 'meter/second'),
    }
