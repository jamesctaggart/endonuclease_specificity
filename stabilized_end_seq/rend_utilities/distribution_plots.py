from plastid.plotting.plotutils import get_fig_axes
from plastid.plotting.plots import get_color_cycle

import copy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import seaborn as sns
import math


def geometric_mean(iterable):
    a = np.array(iterable)
    return a.prod() ** (1.0 / len(a))


def geometric_median(iterable):
    if len(iterable) % 2 == 1:  # Return median value if odd number of entries
        return np.median(iterable)
    else:  # Return geometric mean of middle two values if even number of entries
        return geometric_mean([iterable[(len(iterable) / 2) - 1], iterable[len(iterable) / 2]])


def cdf_plot(data, axes=None, normalizex=False, flip_axes=False, **kwargs):
    """

    :param data: 1D array with data
    :param axes: axis to plot on
    :param logx: if true, plot x on logarithmic scale
    :param xthreshold: factor by which the shaded region will be drawn. Specify None to remove rectangle.
    :param ythreshold: lower and upper limits for highlighting percentiles in data.
    :param normalizex: Normalize CDF position to its median value, allowing for easy superposition of curves
    :return:
    """

    if len(data) == 0:
        return None

    fig, ax = get_fig_axes(axes)

    X = sorted(data)

    if normalizex:
        X = X/geometric_median(X)


    Y = np.arange(0, len(data) + 1) / float(len(data))
    X = np.insert(X, 0, 0.00001)  # This point exists so that there is a line connecting the x axis and your first point

    print len(Y), len(X)
    if 'xlim' not in kwargs:
        xmax = max(X)*2
        xmin = min(X)*0.5
    else:
        xmax = kwargs['xlim'][1]
        xmin = kwargs['xlim'][0]

    if flip_axes:
        ax.step(Y, X, where='post',**kwargs)
    else:
        ax.step(X, Y, where='post',**kwargs)


    return fig, ax