"""
Defines some constants for package

"""

from metpy.units import units

# Stefan-Boltzmann constant in W/m**2/K**4
SIGMA   = 5.670374419e-8

METHODS = [
    'bernard',
    'dimiceli',
    'dimiceli_nws',
    'liljegren',
]

MIN_SPEED          = units.Quantity(2.0, 'knots')
DIMICELI_MIN_SPEED = units.Quantity(1690.0, 'meter per hour')
