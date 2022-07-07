from numpy import exp

def vaporPressure( T ):

  return 6.112 * exp( 17.67 * T / (T + 243.5) )
 
def relativeHumidity( Ta, Td ):

  return vaporPressure( Td ) / vaporPressure( Ta ) * 100.0

def wetBulb( Ta, Td ): 

  RH = relativeHumidity( Ta, Td )

  return -5.806    + 0.672   *Ta -  0.006   *Ta**2       +\
       (  0.061    + 0.004   *Ta + 99.000e-6*Ta**2) * RH +\
       (-33.000e-6 - 5.000e-6*Ta -  1.000e-7*Ta**2) * RH**2

def wetBulbRH( Ta, RH ): 

  return -5.806    + 0.672   *Ta -  0.006   *Ta**2       +\
       (  0.061    + 0.004   *Ta + 99.000e-6*Ta**2) * RH +\
       (-33.000e-6 - 5.000e-6*Ta -  1.000e-7*Ta**2) * RH**2
    
