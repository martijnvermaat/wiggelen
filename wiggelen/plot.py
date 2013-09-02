"""
Plot wiggle tracks.

.. note:: This module depends on the :mod:`matplotlib` package.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import division

import collections
import itertools
import math

from matplotlib import pyplot

from .wiggle import fill


def plot(walker, regions=None, order_by='region', average_threshold=None,
         columns=None, line=True):
    """
    Visualize a wiggle track in a plot.

    For every region, a separate subplot is created.

    :arg walker: Generator yielding tuples of (region, position, value) per
        defined position.
    :type walker: generator(str, int, _)
    :arg regions: Dictionary with regions as keys and (start, stop) tuples as
        values. If not `None`, plot positions from start to stop (both
        including) in these regions. If `None`, plot positions in all
        regions between their first and last defined positions.
    :type regions: dict(str, (int, int))
    :arg order_by: Order the subplots by this key. Possible values are:
        - ``region``: Region name.
        - ``average``: Average value over all defined positions.
        - ``original``: Original position in `walker`.
    :type order_by: str
    :arg average_threshold: Include a subplot for every region with average
        value above or equal to this threshold. If `None`, only include those
        regions with data.
    :type average_threshold: float
    :arg columns: Number of columns to use in subplot arrangement. If `None`,
        a suitable number of columns is chosen automatically.
    :type columns: int
    :arg line: Connect values to create a lineplot. Undefined positions are
        assumed to be `0`.
    :type line: bool

    :return: A tuple containing a matplotlib figure object, a list of
        subplots, the number of rows and the number of columns.
    :rtype: matplotlib.figure.Figure, list(matplotlib.axes.AxesSubplot), int,
        int

    .. note:: This function loads the walker in memory and as such is not
        appropriate for use on very large datasets.
    """
    fig = pyplot.figure(tight_layout=True)
    subplots = collections.OrderedDict()

    if line:
        walker = fill(walker, regions=regions, filler=0, only_edges=True)

    for region, w in itertools.groupby(walker, lambda (r, p, v): r):
        if regions is not None:
            try:
                start, stop = regions[region]
                positions, values = zip(*[(p, v) for r, p, v in w
                                          if start < p <= stop])
            except KeyError:
                continue
        else:
            positions, values = zip(*[(p, v) for r, p, v in w])
            start, stop = min(positions), max(positions)

        average = sum(v for v in values if v is not None) / (stop - start + 1)
        if average_threshold is not None and average < average_threshold:
            continue

        # Temporarily add the subplot at location 111, we'll change that
        # once we know how many subplots there are.
        ax = fig.add_subplot(111, label=region)
        if line:
            ax.plot(positions, values, color='b')
        else:
            ax.plot(positions, values, color='b', linestyle='None', marker=',')
        ax.plot([start, stop], [average, average], color='r', linestyle='--')
        ax.set_xlim(start, stop)
        ax.set_title(region)
        ax.legend([], title='avg=%.2f' % average, loc='best')
        subplots[region] = average, ax

    # If the threshold is zero, also include every region from the defined
    # region that has no data.
    if average_threshold == 0:
        for region in regions:
            if region not in subplots:
                ax = fig.add_subplot(111, label=region)
                ax.set_xlim(*regions[region])
                ax.set_title(region)
                ax.text(0.5, 0.5, 'no data', ha='center', va='center',
                        transform=ax.transAxes, size=24, alpha=0.5)
                subplots[region] = 0, ax

    # Get all the axes in the specified order.
    if order_by == 'region':
        # Alphabetically by region name.
        axes = [subplots[r][1] for r in sorted(subplots)]
    elif order_by == 'average':
        # By average value.
        axes = [subplots[r][1] for r in sorted(subplots, key=subplots.get,
                                               reverse=True)]
    else:
        # By original order in de walker.
        axes = [s[1] for s in subplots.values()]

    columns = columns or int(math.sqrt(len(axes)))
    rows = int((len(axes) - 1) / columns + 1)

    for i, ax in enumerate(axes):
        ax.change_geometry(rows, columns, i + 1)

    return fig, axes, rows, columns
