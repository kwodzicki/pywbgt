from cython.parallel import prange
import numpy

#from libc.stdio cimport printf

cimport numpy
cimport cython

from cliljegren cimport h_sphere_in_air, calc_solar_parameters, calc_wbgt
 
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def chtc( float [:] Tair, float [:] Pair, float [:] speed, float diameter=0.0508  ):
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

  cdef:
    Py_ssize_t i, size = Tair.shape[0]
    float [:] h = numpy.empty( size, dtype = numpy.float32 )                           # Initialize array to write data to

  for i in prange( size, nogil=True ):                                          # Iterate over all values in parallel
    h[i] = h_sphere_in_air( diameter, Tair[i], Pair[i], speed[i] )

  return numpy.asarray(h)                                                     # Reshape to same shape as Tair 

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def test( float [:] Tair, float [:] Pair, float [:] speed, float diameter=0.0508  ):
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
  cdef float [:] hView = h # Initialize array to write data to

  for i in prange( size, nogil=True ):                                          # Iterate over all values in parallel
    hView[i] = h_sphere_in_air( diameter, Tair[i], Pair[i], speed[i] )

  return h                                                     # Reshape to same shape as Tair 


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def solar_parameters( int [::1] year, int [::1] month, int [::1] day, 
                      int [::1] hour, int[::1] minute, int [::1] second,  
                      float [::1] lat, float [::1] lon, solar ):
  """
  Calculate solar parameters based on date and location

  The following quantities are calculated:
    - Modifies solar values to be consistent with normsolar from C function
    - cosine of the solar zenith angle
    - fraction of the solar irradiance due to the direct beam

  Arguments:
    year (ndarray) : Year to compute parameters for
    month (ndarray) : Month to compute parameters for
    day (ndarray) : Day of month
    hour (ndarray) : Hour of day
    minute (ndarray) : Minute of day
    lat (ndarray) : Latitude of location(s) to compute parameters for; decimal
    lon (ndarray) : Longitude of location(s) to compute parameters for; decimal

  Returns:
    tuple : Three (3) ndarrays containing:
      - Potentially modified Solar values
      - cosine of zenith angle
      - fraction of solar irradiance due to the direct beam
 
  """

  cdef:
    int res
    Py_ssize_t i, size = year.shape[0]
    double dday

  out  = solar.reshape( (1, size) ).repeat(3, axis=0)

  cdef float [:, ::1] outView = out

  for i in prange( size, nogil=True ):
    dday  = day[i] + (hour[i]*3600 + minute[i]*60 + second[i])/86400.0          # Compute fractional day
    res   = calc_solar_parameters( year[i], month[i], dday, lat[i], lon[i], 
              &outView[0,i], &outView[1,i], &outView[2,i] )                                             # Run the C function

  return out

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def liljegren( float [:] lat, float [:] lon, int [:] urban,  
               int [:] year, int [:] month, int [:] day, int [:] hour, int [:] minute, int [:] gmt, int[:] avg, 
               float [:] solar, float [:] pres, float [:] Tair, float [:] relhum,
               float [:] speed, float [:] zspeed, float [:] dT): 

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
    lat (ndarray) : Latitude corresponding to data values (decimal)
    lon (ndarray) : Longitude correspondning to data values (decimal)
    urban (ndarray) : Boolean flag indicating if "urban" (1) or
      "rural" (0) for wind speed power law exponent
    year (ndarray) : Year of the data values
    month (ndarray) : Month of the data values
    day (ndarray) : Day of the data values
    hour (ndarray) : Hour of the data values; can be any time zone as long
      as the 'gmt' variable is set appropriately 
    minute (ndarray) : Minute of the data values
    gmt (ndarray) LST-GMT difference  (hours; neative in USA)
    avg (ndarray) : averaging time of the meteorological inputs (minutes)
    solar (ndarray) : solar irradiance (W/m**2)
    pres (ndaray) : barometric pressure (hPa)
    Tair (ndarray) : air (dry bulb) temperature (degree Celsius)
    relhum (ndarray) : Relative humidity (%)
    speed (ndarray) : wind speed (m/s)
    zspeed (ndarray) : height of wind speed measurement (m)
    dT (ndarray) : Vertical temperature difference; upper minus lower
      (degree Celsius)

  Keyword arguments:
    None.

  Returns:
    tuple : The following arrays are returned:
      - Results from the C code sepcifying if there were any errors
      - Estimated 2m wind speed; will be same as input if already 2m temp1
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
    Py_ssize_t i, size = year.shape[0]                                                     # Define size of output arrays based on size of input
    int res
    float est_speed, Tg, Tnwb, Tpsy, Twbg

  out = numpy.empty( (5, size), dtype = numpy.float32 )

  cdef float [:,::1] outView = out 

  # Iterate (in parallel) over all values in the input arrays
  for i in prange( size, nogil=True ):
    res = calc_wbgt( year[i], month[i], day[i], hour[i], minute[i], gmt[i], avg[i],
                               lat[i], lon[i], solar[i], pres[i], Tair[i], relhum[i], speed[i], zspeed[i], dT[i],
                               urban[i], &outView[0,i], &outView[1,i], &outView[2,i], &outView[3,i], &outView[4,i])

  return out
