"""
Package for estimating wetbulb globe temperature

Various algorithms for estimating wetbulb globe
temperature from standard meteorological variables.
 
"""

from .liljegren import wetbulb_globe as liljegrenWBGT
from .bernard   import wetbulb_globe as bernardWBGT
from .dimiceli  import wetbulb_globe as dimiceliWBGT

def wbgt( method, *args, **kwargs ):
    """
    Estimate wet bulb globe temperature

    Wrapper for the various WBGT algorithms provided
    by this package. Set the method to use for 
    estimating WBGT and you're off

    Arguments:
        method (str) : name of the method to use
        *args : Arguments required for each method

    Keyword arguments:
        **kwargs : Passed directly to various algorithms

    Returns:
        dict : Information about wetbulb globe temperature

    """

    method = method.lower()
    if method == 'liljegren':
        return liljegrenWBGT( *args, **kwargs )
    if method == 'bernard':
        return bernardWBGT( *args, **kwargs )
    if method == 'dimiceli':
        return dimiceliWBGT( *args, **kwargs )

    raise Exception( f'Unsupported WBGT method : {method}' )
