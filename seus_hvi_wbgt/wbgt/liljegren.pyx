from cython.parallel import prange
import numpy

#from libc.stdio cimport printf

cimport numpy
cimport cython

from metpy.units import units

from cliljegren cimport *#h_sphere_in_air, calc_solar_parameters, calc_wbgt, Tglobe

from .utils import relative_humidity as rhTd
 
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def chtc( TairIn, PairIn, speedIn, float diameter=0.0508  ):
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

  cdef Py_ssize_t i, size = TairIn.shape[0]
  h = numpy.empty( size, dtype = numpy.float32 )
  cdef:
    float [:] hView = h # Initialize array to write data to
    float [:] Tair  = TairIn.astype( numpy.float32 ) 
    float [:] Pair  = PairIn.astype( numpy.float32 )
    float [:] speed = speedIn.astype( numpy.float32 )

  for i in prange( size, nogil=True ):                                          # Iterate over all values in parallel
    hView[i] = h_sphere_in_air( diameter, Tair[i], Pair[i], speed[i] )

  return h                                                     # Reshape to same shape as Tair 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def solar_parameters( latIn, lonIn, 
                      yearIn, monthIn, dayIn, hourIn, minuteIn, 
                      solar ):
  """
  Calculate solar parameters based on date and location

  The following quantities are calculated:
    - Modifies solar values to be consistent with normsolar from C function
    - cosine of the solar zenith angle
    - fraction of the solar irradiance due to the direct beam

  Arguments:
    lat (ndarray) : Latitude of location(s) to compute parameters for; decimal
    lon (ndarray) : Longitude of location(s) to compute parameters for; decimal
    year (ndarray) : Year to compute parameters for; GMT
    month (ndarray) : Month to compute parameters for; GMT
    day (ndarray) : Day of month; GMT
    hour (ndarray) : Hour of day; GMT
    minute (ndarray) : Minute of day; GMT

  Returns:
    tuple : Three (3) ndarrays containing:
      - Potentially modified Solar values
      - cosine of zenith angle
      - fraction of solar irradiance due to the direct beam
 
  """

  cdef:
    int res
    Py_ssize_t i, size = yearIn.shape[0]
    double dday

  if len( latIn ) == 1:                                                         # If input latitude is only one (1) element, assume lon and urban are also one (1) element and expand all to match size of data
    latIn   = latIn.repeat( size )
    lonIn   = lonIn.repeat( size )

  out  = solar.astype( numpy.float32 ).reshape( (1, size) ).repeat(3, axis=0)

  cdef:
    float [:, ::1] outView = out

    float [::1] lat    = latIn.astype(    numpy.float32 )
    float [::1] lon    = lonIn.astype(    numpy.float32 )
    int   [::1] year   = yearIn.astype(   numpy.int32 ) 
    int   [::1] month  = monthIn.astype(  numpy.int32 )
    int   [::1] day    = dayIn.astype(    numpy.int32 )
    int   [::1] hour   = hourIn.astype(   numpy.int32 )
    int   [::1] minute = minuteIn.astype( numpy.int32 )

  for i in prange( size, nogil=True ):
  #for i in range( size ):
    dday  = day[i] + (hour[i]*60 + minute[i])/1440.0          # Compute fractional day
    res   = calc_solar_parameters( year[i], month[i], dday, lat[i], lon[i], 
              &outView[0,i], &outView[1,i], &outView[2,i] )                                             # Run the C function

  return out

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def globeTemperature( TairIn, TdewIn, PairIn, speedIn, solarIn, fdirIn, czaIn ):
  """
  Compute globe temperature using Liljegren method

  Arguments:
    TairIn (ndarray) : Air (dry bulb) temperature; degree Celsius
    TdewIn (ndarray) : Dew point temperature; degree Celsius
    PairIn (ndarray) : Barometric pressure; hPa
    speedIn (ndarray) : wind speed, m/s
    solarIn (ndarray) : Solar irradiance, W/m**2
    fdirIn (ndarray) : Fraction of solar irradiance due to direct beam
    czaIn (ndarray) : Cosine of solar zenith angle

  Returns:
    ndarray : Globe temperature

  """

  cdef:
    float [:] Tair   = (TairIn + 273.15).astype(  numpy.float32 )
    float [:] Pair   = PairIn.astype(  numpy.float32 )
    float [:] speed  = speedIn.astype( numpy.float32 )
    float [:] solar  = solarIn.astype( numpy.float32 )
    float [:] fdir   = fdirIn.astype(  numpy.float32 )
    float [:] cza    = czaIn.astype(   numpy.float32 )
    float [:] relhum = rhTd(TairIn, TdewIn).astype( numpy.float32 )

    Py_ssize_t i, size = Tair.shape[0]

  out  = numpy.empty( size, dtype=numpy.float32 )

  cdef:
    float [:] outView = out

  for i in prange( size, nogil=True ):
    outView[i] = Tglobe( Tair[i], relhum[i], Pair[i], speed[i], solar[i], fdir[i], cza[i] )

  return out 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.

def psychrometricWetBulb( TairIn, TdewIn, PairIn ):
  """
  Compute psychrometeric wet bulb temperature using Liljegren method

  Arguments:
    TairIn (ndarray) : Air (dry bulb) temperature; degree Celsius
    TdewIn (Quantity) : Dew point temperature; units of temperature
    PairIn (ndarray) : Barometric pressure; hPa

  Returns:
    ndarray : pyschrometric Wet bulb temperature

  """

  cdef:
    float [:] Tair   = (TairIn + 273.15).astype(  numpy.float32 )
    float [:] Pair   = PairIn.astype(  numpy.float32 )
    float [:] relhum = rhTd( TairIn, TdewIn ).astype( numpy.float32 )

    float tmp
    float fill = 0.0
    int   rad  = 0

    Py_ssize_t i, size = Tair.shape[0]

  out  = numpy.full( size, numpy.nan, dtype=numpy.float32 )

  cdef:
    float [:] outView = out

  for i in prange( size, nogil=True ):
    tmp = Twb( Tair[i], relhum[i], Pair[i], fill, fill, fill, fill, rad)
    if tmp > -9999.0: outView[i] = tmp
  return out 


def naturalWetBulb( TairIn, TdewIn, PairIn, speedIn, solarIn, fdirIn, czaIn ):
  """
  Compute natural wet bulb temperature using Liljegren method

  Arguments:
    TairIn (ndarray) : Air (dry bulb) temperature; degree Celsius
    TdewIn (Quantity) : Dew point temperature; units of temperature
    PairIn (ndarray) : Barometric pressure; hPa
    speedIn (ndarray) : wind speed, m/s
    solarIn (ndarray) : Solar irradiance, W/m**2
    fdirIn (ndarray) : Fraction of solar irradiance due to direct beam
    czaIn (ndarray) : Cosine of solar zenith angle

  Returns:
    ndarray : natural Wet bulb temperature

  """

  cdef:
    float [:] Tair   = (TairIn + 273.15).astype(  numpy.float32 )
    float [:] Pair   = PairIn.astype(  numpy.float32 )
    float [:] speed  = speedIn.astype( numpy.float32 )
    float [:] solar  = solarIn.astype( numpy.float32 )
    float [:] fdir   = fdirIn.astype(  numpy.float32 )
    float [:] cza    = czaIn.astype(   numpy.float32 )
    float [:] relhum = rhTd( TairIn, TdewIn ).astype( numpy.float32 )

    float tmp
    int rad = 1

    Py_ssize_t i, size = Tair.shape[0]

  out  = numpy.empty( size, dtype=numpy.float32 )

  cdef:
    float [:] outView = out

  for i in prange( size, nogil=True ):
    tmp = Twb( Tair[i], relhum[i], Pair[i], speed[i], solar[i], fdir[i], cza[i], rad)
    if tmp > -9999.0: outView[i] = tmp

  return out 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def liljegren( latIn, lonIn,
               yearIn, monthIn, dayIn, hourIn, minuteIn,
               solarIn, presIn, TairIn, TdewIn, speedIn,
               urbanIn=None, gmtIn=None, avgIn=None,
               zspeedIn=None, dTIn=None, **kwargs ): 

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
    latIn (ndarray) : Latitude corresponding to data values (decimal).
      Can be one (1) element array; will be expanded to match dates/data
    lonIn (ndarray) : Longitude correspondning to data values (decimal).
      Can be one (1) element array; will be expanded to match dates/data
    yearIn (ndarray) : Year of the data values
    monthIn (ndarray) : Month of the data values
    dayIn (ndarray) : Day of the data values
    hourIn (ndarray) : Hour of the data values; can be any time zone as long
      as the 'gmt' variable is set appropriately 
    minuteIn (ndarray) : Minute of the data values
    solarIn (Quantity) : solar irradiance; units of any power over area
    presIn (Qantity) : barometric pressure; units of pressure
    TairIn (Quantity) : air (dry bulb) temperature; units of temperature
    TdewIn (Quantity) : Dew point temperature; units of temperature
    speedIn (Quatity) : wind speed; units of speed

  Keyword arguments:
    urbanIn (ndarray) : Boolean flag indicating if "urban" (1) or
      "rural" (0) for wind speed power law exponent
      Can be one (1) element array; will be expanded to match dates/data
    gmtIn (ndarray) LST-GMT difference  (hours; neative in USA)
    avgIn (ndarray) : averaging time of the meteorological inputs (minutes)
    zspeedIn (Quantity) : height of wind speed measurement; unit of distance
    dTIn (Quantity) : Vertical temperature difference; upper minus lower;
      unit of temperature

  Returns:
    tuple : The following arrays are returned:
      - Estimated 2m wind speed in meters/second; will be same as input if already 2m temp1
      - Globe temperatures in ndarray
      - Natural wet bulb temperatures in ndarray
      - psychrometric wet bulb temperatures in ndarray 
      - Wet bulb-globe temperatures in ndarray

  Reference: 
    Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp, and R. Sharp:
      Modeling the Wet Bulb Globe Temperature Using Standard Meteorological
      Measurements. The Journal of Occupational and Environmental Hygiene,
      vol. 5:10, pp. 645-655, 2008.

  """

  cdef:
    Py_ssize_t i, size = yearIn.shape[0]                                                     # Define size of output arrays based on size of input
    int res
    float est_speed, Tg, Tnwb, Tpsy, Twbg

  if urbanIn  is None:
    urbanIn  = numpy.zeros( size,       dtype = numpy.int32 )
  if gmtIn    is None:
    gmtIn    = numpy.zeros( size,       dtype = numpy.int32 )
  if avgIn    is None:
    avgIn    = numpy.ones(  size,       dtype = numpy.int32 )
  if zspeedIn is None:
    zspeedIn = numpy.full(  size,  2.0, dtype = numpy.float32 ) * units.meter
  if dTIn     is None:
    dTIn     = numpy.full(  size, -1.0, dtype = numpy.float32 ) * units.degree_Celsius

  if len( latIn ) == 1:                                                         # If input latitude is only one (1) element, assume lon and urban are also one (1) element and expand all to match size of data
    latIn   = latIn.repeat(   size )
    lonIn   = lonIn.repeat(   size )

  if len( urbanIn ) == 1:
    urbanIn = urbanIn.repeat( size )

  cdef:
    float [:] lat   = latIn.astype(   numpy.float32 )
    float [:] lon   = lonIn.astype(   numpy.float32 )
    int   [:] urban = urbanIn.astype( numpy.int32   )

    int [:] year     = yearIn.astype(   numpy.int32 )
    int [:] month    = monthIn.astype(  numpy.int32 )
    int [:] day      = dayIn.astype(    numpy.int32 )  
    int [:] hour     = hourIn.astype(   numpy.int32 )
    int [:] minute   = minuteIn.astype( numpy.int32 )
    int [:] gmt      = gmtIn.astype(    numpy.int32 )
    int [:] avg      = avgIn.astype(    numpy.int32 )

    float [:] solar  =                solarIn.to( 'watt/m**2'      ).magnitude.astype( numpy.float32 )
    float [:] pres   =                 presIn.to( 'hPa'            ).magnitude.astype( numpy.float32 )
    float [:] Tair   =                 TairIn.to( 'degree_Celsius' ).magnitude.astype( numpy.float32 )
    float [:] speed  =                speedIn.to( 'm/s'            ).magnitude.astype( numpy.float32 )
    float [:] zspeed =               zspeedIn.to( 'meter'          ).magnitude.astype( numpy.float32 )
    float [:] dT     =                   dTIn.to( 'degree_Celsius' ).magnitude.astype( numpy.float32 )
    float [:] relhum = (
        100.0 * rhTd( 
            TairIn.to('degree_Celsius').magnitude, 
            TdewIn.to('degree_Celsius').magnitude
        )
    ).astype( numpy.float32 )

  out = numpy.full( (5, size), numpy.nan, dtype = numpy.float32 )

  cdef float [:,::1] outView = out 

  # Iterate (in parallel) over all values in the input arrays
  for i in prange( size, nogil=True ):
    est_speed = speed[i]
    Tg        = 0
    Tnwb      = 0
    Tpsy      = 0
    Twbg      = 0
    res = calc_wbgt( year[i], month[i], day[i], hour[i], minute[i], gmt[i], avg[i],
                     lat[i], lon[i], solar[i], pres[i], Tair[i], relhum[i], speed[i], zspeed[i], dT[i],
                     urban[i], &est_speed, &Tg, &Tnwb, &Tpsy, &Twbg)

    if res == 0:
      outView[0,i] = est_speed
      outView[1,i] = Tg
      outView[2,i] = Tnwb
      outView[3,i] = Tpsy
      outView[4,i] = Twbg

  return {
    'speed' : out[0,:],
    'Tg'    : out[1,:],
    'Tnwb'  : out[2,:],
    'Tpsy'  : out[3,:],
    'Twbg'  : out[4,:]
  }
