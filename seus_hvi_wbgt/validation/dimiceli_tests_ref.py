#!/usr/bin/env python3

import numpy
from pandas import date_range
from metpy.calc import wet_bulb_temperature
from metpy.units import units
from seus_hvi_wbgt.wbgt.dimiceli import *
from seus_hvi_wbgt.wbgt.utils import *


datetime = date_range(
  '2010-09-09T12:00:00', '2010-09-11T12:00:00', freq='1D'
)

lat    = numpy.zeros( datetime.size )
lon    = numpy.zeros( dateime.size )

chfc   = numpy.asarray( [0.315, 0.315, 0.1] )

f_db   = numpy.asarray( [0.67, 0.75, 0.75 ] )
f_dif  = 1.0 - f_db

Tair   = numpy.asarray( [ 30.00,  33.89,  34.44] ) * units( 'degree_Celsius' )
Tdew   = numpy.asarray( [ 20.56,  24.44,  21.11] ) * units( 'degree_Celsius' )
pres   = numpy.asarray( [ 30.08,  29.75,  30.05] ) * units( 'inHg' )
solar  = numpy.asarray( [336.00, 754.00, 579.00] ) * units( 'watt per meter**2')
z      = numpy.asarray( [ 38.44,  36.65,  41.41] ) * units( 'degree' )
speed  = numpy.asarray( [  6.00,   7.00,   3.70] ) * units( 'mile per hour' )
RH     = numpy.asarray( [ 67.50,  54.27,  52.00] )

# These are values calculated by algorithm
Tw_calc   = numpy.asarray( [26.77550, 26.575561, 26.34527] )
Ea_calc   = numpy.asarray( [22.64868, 28.487086, 22.48619] )
Epsa_calc = numpy.asarray( [ 0.897936, 0.9278434, 0.897013] )
B_calc    = numpy.asarray( [ 3614652151, 7098450550, 56176782] )
C_calc	  = numpy.asarray( [ 1197110170, 1309071170, 90440592] )
Tg_calc   = numpy.asarray( [ 33.02, 39.31, 40.65] )
E_calc	  = numpy.asarray( [ 28.54383, 28.584421, 28.24707] )
WBGT_calc = numpy.asarray( [ 28.35, 29.85, 30.02] )

# These are measured truth values
Tg_ref   = numpy.asarray( [ 32.78, 39.44, 40.56] )
WBGT_ref = numpy.asarray( [ 28.30, 29.88, 30.00] )

def relDiff( true, exp ):

  return (exp-true)/true * 100.0

def printError( key, true, exp ):

  print( key )
  for t, e in zip( true, exp ):  
    print( f"  {t:>20.4f} vs. {e:>20.4f}  ---> {relDiff(t, e):10.4f}%" )


cosz  = numpy.cos( z.to('radians').m )

Ea    = atmospheric_vapor_pressure( Tair.to('degree_Celsius').m, 
                        Tdew.to('degree_Celsius').m, 
                        pres.to('hPa'           ).m
)
Epsa  = thermalEmissivity( Tair.to('degree_Celsius').m, 
                           Tdew.to('degree_Celsius').m, 
                           pres.to('hPa'           ).m 
)

B    = factorB( Tair.to('degree_Celsius').m, 
                Tdew.to('degree_Celsius').m, 
                pres.to('hPa'           ).m, 
                solar.to('watt/m**2'    ).m, 
                f_db, f_dif, cosz
)

C    = factorC( speed.to('meter/hour').m, CHFC = chfc ) 

res = dimiceli(lat, lon, datetime,
        solar, pres, Tair, Tdew, speed,
        f_db = f_db,
        cosz = cosz
)

printError( 'Ea',    Ea_calc, Ea )
printError( 'Eps',   Epsa_calc, Epsa )
printError( 'B',     B_calc,    B )
printError( 'C',     C_calc,    C )
printError( 'Globe', Tg_calc, res['Tg'] )
printError( 'Twb',   Tw_calc, res['Tnwb'] )
printError( 'WBGT',  WBGT_calc, res['Twbg'] )
