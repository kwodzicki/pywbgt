from cython.parallel import prange
import numpy

from libc.stdio cimport printf

cimport numpy
cimport cython

from cliljegren cimport *

from metpy.units import units
from pandas import to_datetime, to_timedelta, DatetimeIndex

from .calc import relative_humidity as rhTd
 
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def chtc( Tair, Pair, speed, float diameter=0.0508  ):
  """
  Compute convective heat transfer coefficient 

  Wrapper for the h_sphere_in_air C function from WBGT v1.1

  Arguments:
    Tair (ndarray) : Ambient air temperature in Kelvin
    Pair (ndarray) : Barometric pressure in hPa/mb
    speed (ndarray) : Air/wind seed in m/s

  Keyword arguments:
    diameter (float) : Diameter of the sphere in meters; default value
      is 0.0508 and is taken from the C code

  Returns:
    ndarray : convective heat transfer coefficients

  """

  cdef Py_ssize_t i, size = Tair.shape[0]

  h = numpy.empty( size, dtype = numpy.float32 )

  cdef:
    float [::1] hView     = h # Initialize array to write data to
    float [::1] TairView  = Tair.astype( numpy.float32 ) 
    float [::1] PairView  = Pair.astype( numpy.float32 )
    float [::1] speedView = speed.astype( numpy.float32 )

  for i in prange( size, nogil=True ):                                          # Iterate over all values in parallel
    hView[i] = h_sphere_in_air( 
      diameter, TairView[i], PairView[i], speedView[i]
    )

  return h                                                     # Reshape to same shape as Tair 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def solar_parameters( lat, lon, datetime,
                      solar, avg=None, gmt=None, use_spa=False, **kwargs ):
  """
  Calculate solar parameters based on date and location

  The following quantities are calculated:
    - Modifies solar values to be consistent with normsolar from C function
    - cosine of the solar zenith angle
    - fraction of the solar irradiance due to the direct beam

  Arguments:
    lat (ndarray) : Latitude of location(s) to compute parameters for; decimal
    lon (ndarray) : Longitude of location(s) to compute parameters for; decimal
    datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data

  Keyword arguments:
    avg (ndarray) : averaging time of the meteorological inputs (minutes)
    gmt (ndarray) : LST-GMT difference  (hours; negative in USA)

  Returns:
    tuple : Three (3) ndarrays containing:
      - Potentially modified Solar values
      - cosine of zenith angle
      - fraction of solar irradiance due to the direct beam
 
  """

  if not isinstance( datetime, DatetimeIndex ):
    raise Exception( "The 'datetime' argument must be a 'pandas.DatetimeIndex' object")

  cdef:
    int res, spa = use_spa
    Py_ssize_t i, size = datetime.shape[0]
    double dday

  if avg is None: avg = 1.0
  dt = to_timedelta( avg/2.0, 'minute')
  if gmt is not None: dt = dt + to_timedelta(gmt, 'hour')

  datetime = datetime - dt

  if lat.size <= 1:                                                         # If input latitude is only one (1) element, assume lon and urban are also one (1) element and expand all to match size of data
    lat = lat.repeat( size )
    lon = lon.repeat( size )

  out  = solar.astype( numpy.float32 ).reshape( (1, size) ).repeat(3, axis=0)

  cdef:
    float [:, ::1] outView = out

    float [::1] latView  = lat.astype(    numpy.float32 )
    float [::1] lonView  = lon.astype(    numpy.float32 )

    int [::1] yearView   = datetime.year.values.astype(   numpy.int32 )
    int [::1] monthView  = datetime.month.values.astype(  numpy.int32 )
    int [::1] dayView    = datetime.day.values.astype(    numpy.int32 )  
    int [::1] hourView   = datetime.hour.values.astype(   numpy.int32 )
    int [::1] minuteView = datetime.minute.values.astype( numpy.int32 )
    int [::1] secondView = datetime.second.values.astype( numpy.int32 )
 
  for i in prange( size, nogil=True ):
    #dday  = day[i] + (hour[i]*60 + minute[i] - 0.5*avg[i])/1440.0            # Compute fractional day
    dday  = (
      dayView[i] + 
      ((hourView[i]*60 + minuteView[i])*60 + secondView[i])/86400.0
    )                                                                         # Compute fractional day
    res = calc_solar_parameters(
      yearView[i], monthView[i], dday, latView[i], lonView[i], 
      spa, &outView[0,i], &outView[1,i], &outView[2,i] 
    )                                                                         # Run the C function

  return out

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def globeTemperature( Tair, Tdew, Pair, speed, solar, fdir, cza ):
  """
  Compute globe temperature using Liljegren method

  Arguments:
    Tair (ndarray) : Air (dry bulb) temperature; degree Celsius
    Tdew (ndarray) : Dew point temperature; degree Celsius
    Pair (ndarray) : Barometric pressure; hPa
    speed (ndarray) : wind speed, m/s
    solar (ndarray) : Solar irradiance, W/m**2
    fdir (ndarray) : Fraction of solar irradiance due to direct beam
    cza (ndarray) : Cosine of solar zenith angle

  Returns:
    ndarray : Globe temperature

  """

  cdef:
    float [::1] TairView   = (Tair + 273.15).astype(  numpy.float32 )
    float [::1] PairView   = Pair.astype(  numpy.float32 )
    float [::1] speedView  = speed.astype( numpy.float32 )
    float [::1] solarView  = solar.astype( numpy.float32 )
    float [::1] fdirView   = fdir.astype(  numpy.float32 )
    float [::1] czaView    = cza.astype(   numpy.float32 )
    float [::1] relhumView = rhTd(Tair, Tdew).astype( numpy.float32 )

    Py_ssize_t i, size = Tair.shape[0]

  out = numpy.empty( size, dtype=numpy.float32 )
  cdef float [:] outView = out

  for i in prange( size, nogil=True ):
    outView[i] = Tglobe(
      TairView[i], relhumView[i], PairView[i], 
      speedView[i], solarView[i], fdirView[i], czaView[i]
    )

  return out 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def psychrometricWetBulb( Tair, Tdew, Pair ):
  """
  Compute psychrometeric wet bulb temperature using Liljegren method

  Arguments:
    Tair (ndarray) : Air (dry bulb) temperature; degree Celsius
    Tdew (Quantity) : Dew point temperature; units of temperature
    Pair (ndarray) : Barometric pressure; hPa

  Returns:
    ndarray : pyschrometric Wet bulb temperature

  """

  cdef:
    float [::1] TairView   = (Tair + 273.15).astype(  numpy.float32 )
    float [::1] PairView   = Pair.astype(  numpy.float32 )
    float [::1] relhumView = rhTd( Tair, Tdew ).astype( numpy.float32 )

    float tmp, fill = 0.0
    int   rad  = 0

    Py_ssize_t i, size = Tair.shape[0]

  out = numpy.full( size, numpy.nan, dtype=numpy.float32 )
  cdef float [::1] outView = out

  for i in prange( size, nogil=True ):
    tmp = Twb(
      TairView[i], relhumView[i], PairView[i], fill, fill, fill, fill, rad
    )
    if tmp > -9999.0: outView[i] = tmp
  return out 


def naturalWetBulb( Tair, Tdew, Pair, speed, solar, fdir, cza ):
  """
  Compute natural wet bulb temperature using Liljegren method

  Arguments:
    Tair (ndarray) : Air (dry bulb) temperature; degree Celsius
    Tdew (Quantity) : Dew point temperature; units of temperature
    Pair (ndarray) : Barometric pressure; hPa
    speed (ndarray) : wind speed, m/s
    solar (ndarray) : Solar irradiance, W/m**2
    fdir (ndarray) : Fraction of solar irradiance due to direct beam
    cza (ndarray) : Cosine of solar zenith angle

  Returns:
    ndarray : natural Wet bulb temperature

  """

  cdef:
    float [::1] TairView   = (Tair + 273.15).astype(  numpy.float32 )
    float [::1] PairView   = Pair.astype(  numpy.float32 )
    float [::1] speedView  = speed.astype( numpy.float32 )
    float [::1] solarView  = solar.astype( numpy.float32 )
    float [::1] fdirView   = fdir.astype(  numpy.float32 )
    float [::1] czaView    = cza.astype(   numpy.float32 )
    float [::1] relhumView = rhTd( Tair, Tdew ).astype( numpy.float32 )

    float tmp
    int rad = 1

    Py_ssize_t i, size = Tair.shape[0]

  out = numpy.full( size, numpy.nan, dtype=numpy.float32 )
  cdef float [:] outView = out

  for i in prange( size, nogil=True ):
    tmp = Twb( 
      TairView[i], relhumView[i], PairView[i], 
      speedView[i], solarView[i], fdirView[i], czaView[i], rad
    )
    if tmp > -9999.0: outView[i] = tmp

  return out 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def liljegren( lat, lon, datetime, 
               solar, pres, Tair, Tdew, speed,
               urban=None, gmt=None, avg=None,
               zspeed=None, dT=None, use_spa=False, 
               d_globe=None, **kwargs ): 

  """
  Calculate the outdoor wet bulb-globe temperature

  Cython wrapper for Liljegren C code for calculating the outdoor
  wet bulb-globe temperature,, which is the weighted sum of the\
  air temperature (dry bulb), the globe temperature, and the 
  natural wet bulb temperature:

      Twbg = 0.1 * Tair + 0.7 * Tnwb + 0.2 * Tg.

  The program predicts Tnwb and Tg using meteorological input data
  then combines the results to produce Twbg.

  Arguments:
    lat (ndarray) : Latitude corresponding to data values (decimal).
      Can be one (1) element array; will be expanded to match dates/data
    lon (ndarray) : Longitude correspondning to data values (decimal).
      Can be one (1) element array; will be expanded to match dates/data
    datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
    solar (Quantity) : solar irradiance; units of any power over area
    pres (Qantity) : barometric pressure; units of pressure
    Tair (Quantity) : air (dry bulb) temperature; units of temperature
    Tdew (Quantity) : Dew point temperature; units of temperature
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
    use_spa (bool) : If set, use the National Renewable Energy 
      Laboratory (NREL) Solar Position Algorithm (SPA) to determine
      sun position. Default is to use the build-it, low precision model.
    d_globe (Quantity) : Diameter of the black globe thermometer
      unit of distance

  Returns:
    dict :
      - Tg : Globe temperatures in ndarray
      - Tpsy : psychrometric wet bulb temperatures in ndarray 
      - Tnwb : Natural wet bulb temperatures in ndarray
      - Twbg : Wet bulb-globe temperatures in ndarray
      - Speed : Estimated 2m wind speed in meters/second; will be same as input if already 2m temp1

  Reference: 
    Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp, and R. Sharp:
      Modeling the Wet Bulb Globe Temperature Using Standard Meteorological
      Measurements. The Journal of Occupational and Environmental Hygiene,
      vol. 5:10, pp. 645-655, 2008.

  """

  if not isinstance( datetime, DatetimeIndex ):
    raise Exception( "The 'datetime' argument must be a 'pandas.DatetimeIndex' object")

  # Define size of output arrays based on size of input
  cdef:
    Py_ssize_t i, size = datetime.shape[0]
    int res, spa=use_spa
    float dGlobe, est_speed, Tg, Tnwb, Tpsy, Twbg

  # Set default data averaging interval and compute time offset so that
  # time is in the middle of the sampling interval
  if avg is None: avg = 1.0  
  dt = to_timedelta( avg/2.0, 'minute')

  # If gmt is NOT None (i.e., it is set), then adjust the time delta
  if gmt is not None: dt = dt + to_timedelta(gmt, 'hour')

  # Adjust time using time delta
  datetime = datetime - dt

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
      units.degree_Celsius * numpy.full(
        size, -1.0, dtype=numpy.float32
      )
    ) 

  # Set default black globe thermometer diameter
  if d_globe is None:
    dGlobe = D_GLOBE                                                          # Use default value from C source code
  else:
    dGlobe = d_globe.to('meter').astype(numpy.float32).magnitude              # Ensure is in units of meter, convert to 32-bit float, and get magnitude
    
  if len( lat ) == 1:                                                         # If input latitude is only one (1) element, assume lon and urban are also one (1) element and expand all to match size of data
    lat = lat.repeat( size )
    lon = lon.repeat( size )

  # Define array views for faster/parallel iteration
  cdef:
    float [::1] latView    = lat.astype(   numpy.float32 )
    float [::1] lonView    = lon.astype(   numpy.float32 )
    int   [::1] urbanView  = urban.astype( numpy.int32   )

    int [::1] yearView     = datetime.year.values.astype(   numpy.int32 )
    int [::1] monthView    = datetime.month.values.astype(  numpy.int32 )
    int [::1] dayView      = datetime.day.values.astype(    numpy.int32 )  
    int [::1] hourView     = datetime.hour.values.astype(   numpy.int32 )
    int [::1] minuteView   = datetime.minute.values.astype( numpy.int32 )
    int [::1] secondView   = datetime.second.values.astype( numpy.int32 )
 
    float [::1] solarView  =  solar.to( 'watt/m**2'      ).magnitude.astype( numpy.float32 )
    float [::1] presView   =   pres.to( 'hPa'            ).magnitude.astype( numpy.float32 )
    float [::1] TairView   =   Tair.to( 'degree_Celsius' ).magnitude.astype( numpy.float32 )
    float [::1] speedView  =  speed.to( 'm/s'            ).magnitude.astype( numpy.float32 )
    float [::1] zspeedView = zspeed.to( 'meter'          ).magnitude.astype( numpy.float32 )
    float [::1] dTView     =     dT.to( 'degree_Celsius' ).magnitude.astype( numpy.float32 )
    float [::1] relhumView = (
      100.0 * 
      rhTd( 
        Tair.to('degree_Celsius').magnitude, 
        Tdew.to('degree_Celsius').magnitude
      )
    ).astype( numpy.float32 )

  # Define output array for storing WBGT status code
  out = numpy.full( (5, size), numpy.nan, dtype = numpy.float32 )
  cdef float [:,::1] outView = out 

  # Iterate (in parallel) over all values in the input arrays
  for i in prange( size, nogil=True ):
    # Ensure that variables are defined for each thread
    est_speed = speedView[i]
    Tg        = 0
    Tnwb      = 0
    Tpsy      = 0
    Twbg      = 0

    # Run WBGT code; note that GMT and AVG input arguments are set to
    # zero (0) because all time adjustment is done above in Cython
    res = calc_wbgt( 
      yearView[i], monthView[i], dayView[i], 
      hourView[i], minuteView[i], secondView[i], 
      0, 0,
      latView[i], lonView[i], 
      solarView[i], presView[i], TairView[i], 
      relhumView[i], speedView[i], zspeedView[i], dTView[i],
      urbanView[i], spa, dGlobe, &est_speed, &Tg, &Tnwb, &Tpsy, &Twbg
    )

    # If WBGT was success, then store variables in the outView
    if res == 0:
      outView[0,i] = est_speed
      outView[1,i] = Tg
      outView[2,i] = Tnwb
      outView[3,i] = Tpsy
      outView[4,i] = Twbg

  # Return dict with unit-aware values
  return {
    'Tg'    : units.Quantity(out[1,:], 'degree_Celsius'),
    'Tpsy'  : units.Quantity(out[3,:], 'degree_Celsius'),
    'Tnwb'  : units.Quantity(out[2,:], 'degree_Celsius'),
    'Twbg'  : units.Quantity(out[4,:], 'degree_Celsius'),
    'speed' : units.Quantity(out[0,:], 'meter/second'),
  }
