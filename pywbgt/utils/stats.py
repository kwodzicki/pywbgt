"""
Utilities for basic statistics

"""

import numpy as np

def minmax( *args ):
    """
    Determine minimum and maximum

    The minimum and maximum across all input
    arguments is computed

    Arguments:
        *args : Any number of argumnets

    Returns:
        tuple : Min and Max of input arguments

    """

    out = [float('infinity'), -float('infinity')]
    for arg in args:
        arg = np.asarray(arg)
        out[0] = np.nanmin( [out[0], np.nanmin(arg)] )
        out[1] = np.nanmax( [out[1], np.nanmax(arg)] )
    return out

def mean_abs_error( xdata, ydata ):
    """
    Compute mean absolute error

    Arguments:
        xdata (numpy.ndarray) : Independent variable
        ydata (numpy.ndarray0 : Dependent variable

    Returns:
       float : The mean absolute error

    """

    return np.mean( np.abs( ydata.squeeze() - xdata.squeeze() ), axis=0 )

def root_mean_square_error( xdata, ydata ):
    """
    Compute root mean square error

    Arguments:
        xdata (numpy.ndarray) : Independent variable
        ydata (numpy.ndarray0 : Dependent variable

    Returns:
       float : The root mean square error

    """

    return np.sqrt(
        np.sum( (xdata.squeeze()-ydata.squeeze())**2 ) / xdata.size
    )
