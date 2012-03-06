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

from .wiggle import walk
from .index import index
from .merge import merge


def normalize(*values):
    """
    Normalize values relative to the smallest value.
    """
    min_value = min(*values)
    return [value / min_value for value in values]


def matrix(size, reflexive=False, symmetric=False):
    """
    Create all coordinates in a square matrix.

    :arg size: Width and height of the matrix.
    :type size: int
    :reflexive: Include coordinates on (x, x) diagonal.
    :type reflexive: bool
    :symmetric: Include coordinates (x, y) where x < y.

    :return: All coordinates in the matrix as tuples.
    :rtype: list(tuple(int, int))

    .. todo:: Implement ``reflexive`` and ``symmetric``.
    """
    return [(i, j) for i in range(1, size) for j in range(0, i)]


def distance(*tracks):
    """
    Calculate the pairwise distances between wiggle tracks.

    :arg tracks: List of wiggle tracks.
    :type walkers: list(file)

    :return: Pairwise distances between ``tracks``.
    :rtype: list(int)

    .. todo:: This is not yet finished. For one thing, we need the summed
        values per track, which we don't store in the index at the moment.
    .. todo:: Noise reduction with threshold.
    .. todo:: Choose pairwise comparison function.
    .. todo:: Can we refactor this by generalizing to a merge.pairwise?
    .. todo:: Not creating the index is not an option here, update note on
        assumption for this.
    .. todo:: Check where this goes wrong if we cannot .seek() the tracks.
    """
    walkers = [walk(track, force_index=True) for track in tracks]

    # We construct a list of comparisons for the merger, where each comparison
    # is a tuple of (left, right, weight_left, weight_right).
    comparisons = []
    coverages = [index(track, force=True)[0]['sum'] for track in tracks]
    for left, right in matrix(len(tracks)):
        weight_left, weight_right = normalize(coverages[left], coverages[right])
        comparisons.append( (left, right, weight_left, weight_right) )

    # Pairwise distance function
    distance_function = lambda x, y: abs(x - y) / ((x + 1) * (y + 1))

    # The merger returns a matrix of values, one for each pairwise comparison.
    def merger(values):
        results = {}
        for left, right, weight_left, weight_right in comparisons:
            value_left, value_right = values[left], values[right]
            if value_left is None and value_right is None:
                result = None
            else:
                x = weight_left * value_left if value_left else 0
                y = weight_right * value_right if value_right else 0
                result = distance_function(x, y)
            results[left, right] = result
        return results

    # Sum the results.
    counters = defaultdict(lambda: 0)
    denominators = defaultdict(lambda: 1)
    for _, _, values in merge(*walkers, merger=merger):
        for comparison, value in values.items():
            if value is not None:
                counters[comparison] += value
                denominators[comparison] += 1

    # Distance matrix.
    distances = {}
    for comparison in counters:
        distances[comparison] = counters[comparison] / denominators[comparison]
    return distances
