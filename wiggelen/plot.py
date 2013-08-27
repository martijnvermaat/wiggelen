"""
Plot wiggle tracks.

.. note:: This module depends on the :mod:`matplotlib` package.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import division

import itertools
import math

from matplotlib import pyplot


def plot(walker, columns=None):
    """
    Visualize a wiggle track in a plot.

    For every region, a separate subplot is created.

    :arg walker: Generator yielding tuples of (region, position, value) per
        defined position.
    :type walker: generator(str, int, _)
    :arg columns: Number of columns to use in subplot arrangement. If `None`,
        a suitable number of columns is chosen automatically.

    :return: A tuple containing a matplotlib figure object, a list of
        subplots, the number of rows and the number of columns.
    :rtype: matplotlib.figure.Figure, list(matplotlib.axes.AxesSubplot), int,
        int
    """
    fig = pyplot.figure(tight_layout=True)
    axes = []

    for region, w in itertools.groupby(walker, lambda (r, p, v): r):
        # Todo: Use wiggelen.fill to have None on undefined positions.
        positions, values = zip(*[(p, v) for r, p, v in w])

        # Temporarily add the subplot at location 111, we'll change that
        # once we know how many subplots there are.
        ax = fig.add_subplot(111, label=region)
        ax.plot(positions, values)
        ax.set_xlim(min(positions), max(positions))
        ax.set_title(region)
        axes.append(ax)

    columns = columns or int(math.sqrt(len(axes)))
    rows = int((len(axes) - 1) / columns + 1)

    for i, ax in enumerate(axes):
        ax.change_geometry(rows, columns, i + 1)

    return fig, axes, rows, columns
