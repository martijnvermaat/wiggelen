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

from collections import defaultdict

from .merge import merge


def getweights(coverage_left, coverage_right):
    min_coverage = min(coverage_right, coverage_left)
    return coverage_left / min_coverage, coverage_right / min_coverage


def pairwise_comparisons(count):
    return [(i, j) for i in range(1, count) for j in range(0, i)]


def distance(*walkers):
    """
    Calculate the pairwise distances between wiggle tracks.

    This assumes the walkers have their regions in the same order. You can
    force this by using indices. Example::

        >>> from wiggelen import walk
        >>> from wiggelen.distance import distance
        >>> walkers = [walk(track, force_index=True) for track in tracks]
        >>> distance(*walkers)

    :arg walkers: List of generators yielding tuples of (region, position,
        value) per defined position.
    :type walkers: list(generator(str, int, _))

    :return: Pairwise distances between ``walkers``.
    :rtype: list(int)

    .. todo:: This is not yet finished. For one thing, we need the summed
        values per track, which we don't store in the index at the moment.
    .. todo:: Noise reduction with threshold.
    .. todo:: Choose pairwise comparison function.
    """
    # Temporary values, eventually this should be stored in the index.
    coverages = [23688, 9105, 34243, 34243, 34243][:len(walkers)]

    comparisons = pairwise_comparisons(len(walkers))

    weights = dict(((left, right), getweights(coverages[left], coverages[right]))
                   for left, right in comparisons)

    def merger(values):
        results = {}
        for left, right in comparisons:
            x = weights[left, right][0] * values[left] if values[left] else 0.0
            y = weights[left, right][1] * values[right] if values[right] else 0.0
            results[left, right] = float(abs(x - y)) / ((x + 1) * (y + 1))
        return results

    counters = defaultdict(lambda : 0)
    denominators = defaultdict(lambda : 1)
    for _, _, values in merge(*walkers, merger=merger):
        for comparison, value in values.items():
            counters[comparison] += value
            denominators[comparison] += 1  # Todo: only count 1 where at least one of the two is defined

    return dict((c, counters[c] / denominators[c]) for c in comparisons)
