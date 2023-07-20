"""
Scatter plot with color binned symbols

Create scatter plot where the symbols are colored based on what
bin of a variable they fall in. For instances, can color based on
wind speed, with speeds of 0-2 m/s red, 2-4 m/s blue, etc.

"""

from pandas import DataFrame
from numpy import random, arange

from matplotlib import pyplot as plt
from matplotlib import colors as mcolors

from ..stats import mean_abs_error
from .utils import bin_var, thin_data, textbf

def plot_by_hue(xdata, ydata, hue,
        bin_size     = None,
        bin_max      = None,
        xlabel       = '',
        ylabel       = '',
        mae_fontsize = 'small',
        ax           = None,
        figsize      = (6,6),
        seed         = 42,
        npoints_max  = 10**5,
        **kwargs,
    ):
    """
    Scatter plot data and color symbols

    Arguments:
        xdata (numpy.ndarray) : Values for the x-axis
        ydata (numpy.ndarray) : Values for the y-axis
        hue (numpy.ndarray) : Data to bin by for coloring

    Keyword arguments:
        bin_size (int,float) : Size of the bins for coloring
        bin_max (int,float) : Right bound of largest bin. Note that
            'bin_min' is always set to zero (0)
        xlabel (str) : Label for the x-axis
        ylabel (str) : Label for the y-axis
        mae_fontsize (str) : Size of the font to use when drawing
            mean absolute error for data.
        ax (matplotlib.pyplot.Axes) : Axis to plot the data on.
        figsize (tuple) : Size of the figure; ignored if ax= is used.
        seed (int) : Seed for random number generator used to supbset
            values for plotting

    Returns:
        tuple : Reference to the figure and axis used in the ploat

    """

    data = DataFrame(
        {
            'x'   : xdata,
            'y'   : ydata,
            'hue' : bin_var(hue, bin_size, bin_max),
        }
    )

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot( 1, 1, 1 )
    else:
        fig = ax.get_figure()

    colors = list( mcolors.TABLEAU_COLORS.keys() )
    colors = {
        key : colors[i % len(colors)]
        for i, key in enumerate(data['hue'].unique())
    }

    for i, (hue_key, hue_df) in enumerate(data.groupby('hue')):
        xvals = hue_df['x']
        yvals = hue_df['y']

        xvals, yvals = thin_data(
            xvals, yvals,
            npoints_total = xdata.size,
            npoints_max   = npoints_max,
            seed          = seed,
        ) 

        ax.scatter(
            xvals, yvals, 
            marker = '.',
            color  = colors[hue_key],
            label  = hue_key,
            **kwargs,
        )
        try:
            mae = mean_abs_error(hue_df['x'], hue_df['y'])
        except:
            pass
        else:
            ax.text(0.98, 0.07*i+0.05, textbf( f"MAE = {mae:0.3f}" ),
                horizontalalignment='right',
                transform = ax.transAxes,
                color     = colors[hue_key],
                weight    = 'bold',
                fontsize  = mae_fontsize,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)


    return fig, ax
