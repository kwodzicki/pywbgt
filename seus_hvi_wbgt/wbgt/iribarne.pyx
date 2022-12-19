from cython.parallel import prange
import numpy

cimport numpy
cimport cython

cdef:
  float C0  =   26.66082
  float C1  =    0.0091379024
  float C2  = 6106.396

  float Cp  = 1004.0    # J/K/kg
  float L   =    2.54e6 # J/kg
  float eps =    0.622
  float f   = Cp/L/eps

cdef extern from "math.h":
    
  double exp(double x) nogil

cdef float abs( float val ) nogil:

  if (val < 0.0):
    return (-1) * val
  return val

cdef float saturation( float T ) nogil:

  return <float>exp( <float>C0 - C1*T - C2/T )

cdef float Tw( float Ta, float Td, float p, int maxfev, float xtol) nogil:
  """
  Iteratively compute the wet-bulb temperature 

  Uses the Iribarne and Godson (1981) method to calculate the
  wet-bulb temperature given dry-bulb temperature, dew point
  temperature, and atmospheric pressure.

  Arguments:
    Ta : Dry-bulb temperature (Kelvin)
    Td : Dew-point temperature (Kelvin)
    p  : Atmoshperic pressure (hPa)
    maxfev : Maximum number of iterations to perform
    xtol : Calculation will terminate if the relative error between
      two consecutive iterates is at most xtol

  Returns:
    float : Wet-bulb temperature in Kelvin

  """

  cdef float es, ed, ew, s, Tw, de, der, dd

  es  = saturation( Ta )
  ed  = saturation( Td )
  s   = (es-ed)/(Ta-Td)

  Tw  = (Ta*f*p + Td*s)/(f*p + s)

  while maxfev > 0:                                                             # While have NOT reached max iteration
    ew   = saturation( Tw )                                                     # Compute saturation vapor pressure for wet-bulb temperature
    de   = f*p*(Ta-Tw) - (ew-ed)                                                # Computer error
    der  = ew*(C1-C2/Tw**2) - f*p
    dd   = de / der
    Tw  -= dd
    if abs(dd) < xtol: return Tw
    maxfev -= 1

  return 0.0/0.0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def psychrometricWetBulb(Ta, Td, p=None, **kwargs):
  """

  Arguments:
    Ta : Dry-bulb temperature (Kelvin)
    Td : Dew-point temperature (Kelvin)

  Keyword arguments:
    p  : Atmoshperic pressure (hPa)
    maxfev : Maximum number of iterations to perform
    xtol : Calculation will terminate if the relative error between
      two consecutive iterates is at most xtol

  """

  cdef:
    Py_ssize_t i, size = Ta.size
    int maxfev = kwargs.get('maxfev', 25)
    float xtol = kwargs.get('xtol',   0.02)

  out = numpy.empty( size, dtype = numpy.float32 )
  cdef:
    float [::1] TaView  = Ta.astype( numpy.float32 )
    float [::1] TdView  = Td.astype( numpy.float32 )
    float [::1] pView   = numpy.full( size, 1013.25, dtype=numpy.float32 ) \
      if p is None else p.astype(  numpy.float32 )
    float [::1] outView = out


  for i in prange( size, nogil=True ):
    outView[i] = Tw( 
      TaView[i]+273.15, TdView[i]+273.15, pView[i], maxfev, xtol
    ) - 273.15

  return out
