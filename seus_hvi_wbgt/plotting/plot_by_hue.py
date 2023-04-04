"""

"""
import warnings

from pandas import DataFrame
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors

from ..stats import mean_abs_error
from .utils import bin_var, flip, textbf

def plot_by_hue(xdata, ydata, hue,
        bin_size     = None,
        bin_max      = None,
        xlabel       = '',
        ylabel       = '',
        title        = '',
        mae_fontsize = 'small', 
        ax           = None,
        figsize      = (6,6),
        **kwargs,
    ):
    """
    bin_label : A string for the column to use to color the points.
        Can be a tuple where zeroth element is column name and first
        element is units. This name will be used for the legend title,
        with units added to name if they are input

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

    for jj, (hue_key, hue_df) in enumerate(data.groupby('hue')):
        ax.scatter(
            hue_df['x'], hue_df['y'],
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
            ax.text(0.98, 0.07*jj+0.05, textbf( f"MAE = {mae:0.3f}" ),
                horizontalalignment='right',
                transform = ax.transAxes,
                color     = colors[hue_key],
                weight    = 'bold',
                fontsize  = mae_fontsize,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)


    return fig, ax

def main( df, x_long_name, x_short_name, y_long_name, y_short_name,
        hues = [ ('Elevation', 'm'), ('Irradiance', 'W/m$^{2}$'), ('Wind', 'mph'), 'Day-Night'],
        methods = ['liljegren', 'dimiceli'],
    ):

    """
    Arguments:
        df (pandas.DataFrame) : All data that will be used for plotting
        var_long_name (str) : Long-name of the variable to plot
        var_short_name (str) : Short name of the variable to plot

    """

    for hue in hues:
        fig, axes = plt.subplots( 1, len(methods), figsize=(6,4) )
        for ax, method in zip(axes, methods):
            title = [y_long_name, f'{method} versus {x_long_name}']
            _ = plotHue(
                subset, var_short_name, f'{var_short_name} {method}',
                ax        = ax,
                xlabel    = f'ECONet {var_short_name} [K]',
                ylabel    = f'{method} {var_short_name} [K]',
                title     = os.linesep.join( title ),
                hue       = hue,
                **figkwargs,
            )

        fig.suptitle(y_long_name)
        fig.subplots_adjust(**pos)
        makeSquare( axes, oneone=True )

        fname = ' '.join( title ) + f'_{hue[0]}.png'
        if saveplot: fig.savefig( os.path.join( figdir, fname ), dpi=dpi )
