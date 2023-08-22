import numpy
from metpy.units import units
from matplotlib import pyplot as plt

HCONST = 5.3865e-8


cosz = numpy.cos( numpy.asarray([38.44, 36.65, 41.41]) * numpy.pi / 180.0  )
cosz = numpy.concatenate( (cosz, [0.9244, 0.9689, 0.9226]) )

S    = numpy.asarray( [336, 754, 579, 448.7, 173, 387.98] ) 
C    = numpy.asarray( [1197110170, 1309071170, 90440592, 444878843.2, 378649115.5, 599274858.9] )
u    = numpy.asarray( [6, 7, 3.7, 2.164, 2.734, 3.86] ) * units('mile/hour')
u    = u.to('meter/hour').m

H    = numpy.asarray( [0.1325695, 0.098524678, 0.127660528] )
h    = HCONST * C / u**0.58

def func( x, a, b, c ):

  return a * x[0]**b * x[1]**c

from scipy.optimize import curve_fit

init_vals = [0.5, 0.5, 0.5]
init_vals = [2, 3, 1.5]
xdata     = [S, cosz]
ydata     = h 
popt, pcov = curve_fit( func, xdata, ydata, p0 = init_vals )

print( popt )
print( pcov )

print( h )
print( func(xdata, *popt) )

print( h[3:] - H )

#fig, ax = plt.subplots(1, 2)
#ax[0].scatter( S, h )
#ax[1].scatter( cosz, h )
#plt.show()
