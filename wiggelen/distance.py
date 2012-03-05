"""
Calculate the distance between two wiggle tracks using a metric designed for
multisets.

.. todo: Perhaps we should just put all this in the wiggelen.merge module.

.. Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
.. Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
.. Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import division

from .merge import merge


def weights(coverage_left, coverage_right):
    if coverage_left < coverage_right:
        return 1.0, coverate_right / coverage_left

    return coverage_left / coverage_right, 1.0


def distance(left, right):
    """
    Calculate the distance between two wiggle tracks.

    This assumes the walkers have their regions in the same order. You can
    force this by using indices. Example::

        >>> from wiggelen import walk
        >>> from wiggelen.distance import distance
        >>> left = walk(left_track, force_index=True)
        >>> right = walk(right_track, force_index=True)
        >>> distance(left, right)

    :arg left: Generator yielding triples of (region, position, value)
        per defined position.
    :type left: generator(str, int, float)
    :arg right: Generator yielding triples of (region, position, value)
        per defined position.
    :type right: generator(str, int, float)

    :return: Distance between ``left`` and ``right``.
    :rtype: int

    .. todo:: This is not yet finished. For one thing, we need the summed
        values per track, which we don't store in the index at the moment.
    .. todo:: Noise reduction with threshold.
    .. todo:: Choose pairwise comparison function.
    """
    coverage_left = 23688
    coverage_right = 9105

    weight_left, weight_right = weights(coverage_left, coverage_right)

    def merger(values):
        x = weight_left * values[0] if values[0] else 0.0
        y = weight_right * values[1] if values[1] else 0.0
        return float(abs(x - y)) / ((x + 1) * (y + 1))

    counter = 0
    denominator = 1

    for _, _, value in merge(left, right, merger=merger):
        counter += value
        denominator += 1

    return counter / denominator
