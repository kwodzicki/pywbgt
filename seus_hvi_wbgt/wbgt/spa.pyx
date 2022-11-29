from cython.parallel import prange
import numpy

from libc.stdio cimport printf

cimport numpy
cimport cython

from cspa cimport *

SPA_ZA     = 0
SPA_ZA_INC = 1
SPA_ZA_RTS = 2
SPA_ALL    = 3

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False)   # Deactivate initialization checking.
def solar_parameters( latIn, lonIn, 
                      yearIn, monthIn, dayIn, hourIn, minuteIn, secondIn, 
  function       = SPA_ALL,
  timezone       = None, 
  delta_ut1      = None,
  delta_t        = None,
  elevation      = None,
  pressure       = None,
  temperature    = None,
  slope          = None,
  azm_rotation   = None,
  atmos_refract  = None,
  ):

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

  Keyword arguments:
    avgIn (ndarray) : averaging time of the meteorological inputs (minutes)

  Returns:
    tuple : Three (3) ndarrays containing:
      - Potentially modified Solar values
      - cosine of zenith angle
      - fraction of solar irradiance due to the direct beam
 
  """

  cdef:
    int res
    Py_ssize_t i, size = yearIn.shape[0]
    spa_data spa

  out   = numpy.empty( (8, size,), dtype=numpy.float64 )
  zeros = numpy.zeros( (size,), dtype=numpy.float64 )

  if timezone      is None: timezone      = zeros 
  if delta_ut1     is None: delta_ut1     = zeros
  if delta_t       is None: delta_t       = zeros
  if elevation     is None: elevation     = zeros
  if pressure      is None: pressure      = zeros
  if temperature   is None: temperature   = zeros
  if slope         is None: slope         = zeros
  if azm_rotation  is None: azm_rotation  = zeros
  if atmos_refract is None: atmos_refract = zeros

  cdef:
    int func = function

    int    [::1] year              = yearIn.astype(   numpy.int32 ) 
    int    [::1] month             = monthIn.astype(  numpy.int32 )
    int    [::1] day               = dayIn.astype(    numpy.int32 )
    int    [::1] hour              = hourIn.astype(   numpy.int32 )
    int    [::1] minute            = minuteIn.astype( numpy.int32 )
    double [::1] second            = secondIn.astype( numpy.float64 )
    double [::1] lat               = latIn.astype(    numpy.float64 )
    double [::1] lon               = lonIn.astype(    numpy.float64 )
    double [::1] timezoneView      = timezone
    double [::1] delta_ut1View     = delta_ut1
    double [::1] delta_tView       = delta_t
    double [::1] elevationView     = elevation
    double [::1] pressureView      = pressure
    double [::1] temperatureView   = temperature 
    double [::1] slopeView         = slope 
    double [::1] azm_rotationView  = azm_rotation
    double [::1] atmos_refractView = atmos_refract

    double [:,::1] outView = out

  for i in prange( size, nogil=True ):
    spa = spa_data()

    spa.function      = func
    spa.year          = year[i]
    spa.month         = month[i]
    spa.day           = day[i]
    spa.hour          = hour[i]
    spa.minute        = minute[i]
    spa.second        = second[i]
    spa.latitude      = lat[i]
    spa.longitude     = lon[i]
    spa.timezone      = timezoneView[i]
    spa.delta_ut1     = delta_ut1View[i]
    spa.delta_t       = delta_tView[i]
    spa.elevation     = elevationView[i]
    spa.pressure      = pressureView[i]
    spa.temperature   = temperatureView[i]
    spa.slope         = slopeView[i]
    spa.azm_rotation  = azm_rotationView[i]
    spa.atmos_refract = atmos_refractView[i]

    #printf( "[%d] timezone : %f\n", i, spa.timezone )
    #printf( "[%d] delta_ut1   : %f\n", i, spa.delta_ut1     )
    #printf( "[%d] delta_t     : %f\n", i, spa.delta_t       )
    #printf( "[%d] elevation   : %f\n", i, spa.elevation     )
    #printf( "[%d] pressure    : %f\n", i, spa.pressure      )
    #printf( "[%d] temperature : %f\n", i, spa.temperature   )
    #printf( "[%d] slope       : %f\n", i, spa.slope         )
    #printf( "[%d] azm_rot     : %f\n", i, spa.azm_rotation  )
    #printf( "[%d] atm_ref     : %f\n", i, spa.atmos_refract )

    res = spa_calculate( &spa )
    if (res == 0):
      outView[0,i] = spa.zenith        # topocentric zenith angle [degrees]
      outView[1,i] = spa.azimuth_astro # topocentric azimuth angle (westward from south) [for astronomers]
      outView[2,i] = spa.azimuth       # topocentric azimuth angle (eastward from north) [for navigators and solar radiation]
      outView[3,i] = spa.incidence     # surface incidence angle [degrees]
      outView[4,i] = spa.suntransit    # local sun transit time (or solar noon) [fractional hour]
      outView[5,i] = spa.sunrise       # local sunrise time (+/- 30 seconds) [fractional hour]
      outView[6,i] = spa.sunset        # local sunset time (+/- 30 seconds) [fractional hour]
      outView[7,i] = spa.r             # earth radius vector [Astronomical Units, AU]
      #printf("Julian Day:    %.6f\n",spa.jd);
      #printf("L:             %.6e degrees\n",spa.l);
      #printf("B:             %.6e degrees\n",spa.b);
      #printf("R:             %.6f AU\n",spa.r);
      #printf("H:             %.6f degrees\n",spa.h);
      #printf("Delta Psi:     %.6e degrees\n",spa.del_psi);
      #printf("Delta Epsilon: %.6e degrees\n",spa.del_epsilon);
      #printf("Epsilon:       %.6f degrees\n",spa.epsilon);
      #printf("Zenith:        %.6f degrees\n",spa.zenith);
      #printf("Azimuth:       %.6f degrees\n",spa.azimuth);
      #printf("Incidence:     %.6f degrees\n",spa.incidence);

  return {
    'zenith'        : out[0,:], 
    'azimuth_astro' : out[1,:], 
    'azimuth'       : out[2,:], 
    'incidence'     : out[3,:], 
    'suntransit'    : out[4,:], 
    'sunrise'       : out[5,:], 
    'sunset'        : out[6,:],
    'r'             : out[7,:],
  } 
