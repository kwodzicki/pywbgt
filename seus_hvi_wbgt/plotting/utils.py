"""
Utilities for making plots pretty

"""

import string
import itertools

import numpy as np

BOTTOM = {'y' : 0.05, 'verticalalignment'   : 'bottom'}
TOP    = {'y' : 0.95, 'verticalalignment'   : 'top'}
LEFT   = {'x' : 0.05, 'horizontalalignment' : 'left'}
RIGHT  = {'x' : 0.95, 'horizontalalignment' : 'right'}

def bin_var( data, bin_size, bin_max=None ):
    """
    Create bin labels for data

    Given an array of values, and bin size, and a
    maximum bin bound, a bin label for each data point is
    created. This is done for plotting purposes so that data
    can be grouped by the bin labels and color-coded.

    Arguments:
        data (numpy.ndarray) : Data to bin
        bin_size (int,float) : The size of the bins
        bin_max (int,float) : The right bound fo the last (largest)
            bin.

    Returns:
        numpy.ndarray : Numpy array containing string labels
            for which bin each data point falls in

    """

    if bin_size is None:
        return data

    if np.any( ~np.isfinite(data) ):
        warnings.warn(
            "Non-finite values in binnig data. "
            "These values will end up in LAST (largest) bin!",
            RuntimeWarning,
        )

    if isinstance(bin_size, (list,tuple,np.ndarray)):
        bins = np.asarray( bin_size )
    elif bin_max is not None:
        bins    = np.arange(0, bin_max+1, bin_size )
    else:
        raise Exception( 'Failed to construct bins for binning!' )

    idx     = np.digitize(data, bins)

    bbins  = [f"{bins[i-1]:4d} - {bins[i]:4d}" for i in range(1, len(bins))]
    bbins.append( f"> {bins.max()}" )

    return np.asarray(bbins)[idx-1]

def flip(items, ncol):
    """
    To flip plot legend

    Arguments:
        items (list) : Items for the legend
        ncol (int) : Number of columns in the legend

    Returns:
        itertools.chain

    """

    return itertools.chain(*[items[i::ncol] for i in range(ncol)])


def textbf( txt ):
    """
    To make text bold face using LaTeX rendering

    Arguments:
        txt (str) : String to make bold face

    Returns:
        str : LaTeX formatted bold face string

    """

    return r"\textbf{"+txt+"}"

def make_square( axes, oneone=False ):
    """
    Make axes square

    Arguments:
        axes (list) : matplotib axes objects to make
            square.

    Keyword arguments:
        oneone (bool) : If set, will draw a one-to-one line
            on every axis

    Returns:
        None

    """

    ax_min, ax_max = float('inf'), -float('inf')
    for axis in axes:
        y_min, y_max = axis.get_ylim()
        x_min, x_max = axis.get_xlim()
        ax_min = min(ax_min, y_min, x_min)
        ax_max = max(ax_max, y_max, x_max)
    lim = [ax_min, ax_max]
    for axis in axes:
        axis.set_ylim( lim )
        axis.set_xlim( lim )
        axis.set_aspect('equal')
        if oneone and axis.axison:
            axis.plot( lim, lim, zorder=0, color='black', linestyle='-' )

def get_all_handles_labels( axes ):
    """
    Get all labels across axes

    When using the plot_by_hue() function, some panels in multipanel
    plot may not have all the labels of other plots. Thus, this
    function will iterate over all the axes in the numpy array and
    determine which axis has the most handles/labels. The axis along
    with the handles and labels for this 'biggest' axis are returned.

    Arguments:
        axes (numpy.ndarray) : Array of axes to iterate over

    Keywords arguments:
        None.

    Returns:
        tuple : The axis that contained the most labels, the
            handles for the legend, and the labels for the legend.

    """

    nlabels =  0
    for axis in axes.flatten():
        handles_tmp, labels_tmp = axis.get_legend_handles_labels()
        if len(handles_tmp) <= nlabels:
            continue
        nlabels = len( handles_tmp )
        handles, labels = handles_tmp, labels_tmp
        axis_all = axis

    return axis_all, handles, labels

def add_legend( source_ax, title='', ncol=5, legend_ax=None ):
    """
    Arguments:
        source_ax (matplotlib.Axes) : Source axis to get the
            data label information from.

    Keyword arguments:
        legend_ax (matplotlib.Axes) : Axis to draw legend in.
            Defaults to figure of source_ax

    """

    if isinstance(source_ax, np.ndarray):
        source_ax, handles, labels = get_all_handles_labels( source_ax )

    if legend_ax is None:
        legend_ax = source_ax.get_figure()
    else:
        legend_ax.set_axis_off()

    lgnd = legend_ax.legend(
        flip(handles, ncol), flip(labels, ncol),
        loc      = 'lower center',
        ncol     = ncol,
        title    = title,
        fontsize = 'medium',
        columnspacing = 0.5,
        handletextpad = 0.1,
    )
    for ll in lgnd.legendHandles:
        ll._sizes=[80]


def axes_labels(axes,
        location  = 'top left',
        formatter = '({})',
        uppercase = False,
        rows      = False
    ):
    """
    Add labels to each axes of multipanel

    Given an array of axes, letter labels will be written
    on plot to distinguish the panels

    Arguments:
        axes (numpy.ndarray) : Array of axes as returned by a call to
            matplotlib.pyplot.subplots() or other method of adding
            axes to a plot.

    Keyword arguments:
        location (str) : Location of the label. Valid options are
            'top left', 'top right', 'bottom left', and 'bottom right'.
        formatter (str) : A formatting string for how label should
            be formatted. Default is '(a)', '(b)', etc.
        uppercase (bool) : If set, use uppercase letters
        rows (bool) : If set, write letters across rows before
            going across columns.

    Returns:
        None.

    """

    loc = BOTTOM if 'bottom' in location else TOP
    loc.update(LEFT if 'left' in location else RIGHT)
    if len(loc) != 4:
        raise Exception( 'Error in location!')
 
    letters  = string.ascii_uppercase if uppercase else string.ascii_lowercase

    if axes.ndim == 2:
        if rows:
            axes = axes.T
        axes = axes.flatten()
    
    for axis, letter in zip(axes, letters):
        axis.text(
            s         = formatter.format(letter), 
            transform = axis.transAxes,
            **loc,
        )
