"""
Calculate the distance between two wiggle tracks using a metric designed for
multisets.

This module implements the algorithm from the `wiggledist <https://humgenprojects.lumc.nl/trac/wiggledist/>`_
program, an efficient tool to assess similarity of next generation sequencing
datasets.

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


# Todo: Give these metrics a name.
_metric_a = lambda x, y : abs(x - y) / ((x + 1) * (y + 1))
_metric_b = lambda x, y : abs(x - y) / (x + y + 1)
_metric_c = lambda x, y : (max(x, y) * abs(x - y)) \
                          / (((x * x) + 1) * ((y * y) + 1))
_metric_d = lambda x, y : abs(x - y) / (max(x, y) + 1)


#: Predefined pairwise distance metrics.
metrics = {'a': _metric_a,
           'b': _metric_b,
           'c': _metric_c,
           'd': _metric_d}


def normalize(*values):
    """
    Normalize values relative to the smallest value.

    :arg values: List of values.
    :type values: list(float)

    :return: Scale the values such that the minimum is 1.
    :rtype: list(float)
    """
    min_value = min(*values)
    if not min_value:
        return [1 for _ in values]
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
    return [(i, j) for i in range(1, size) for j in range(i)]


def distance(*tracks, **options):
    """
    Calculate the pairwise distances between wiggle tracks.

    :arg tracks: List of wiggle tracks.
    :type walkers: list(file)
    :keyword metric: Pairwise distance metric (default: a).
    :type merger: function(float, float -> float)

    :return: Pairwise distances between ``tracks`` as a mapping from
        coordinates in the distance matrix to their values.
    :rtype: dict((int, int), float)

    .. todo:: This is not yet finished. For one thing, we need the summed
        values per track, which we don't store in the index at the moment.
    .. todo:: Noise reduction with threshold.
    .. todo:: Choose pairwise comparison function.
    .. todo:: Can we refactor this by generalizing to a merge.pairwise?
    .. todo:: Not creating the index is not an option here, update note on
        assumption for this.
    .. todo:: Check where this goes wrong if we cannot .seek() the tracks.
    """
    metric = options.get('metric', metrics['a'])

    # We construct a list of comparisons for the merger, where each comparison
    # is a tuple of (left, right, weight_left, weight_right).
    comparisons = []
    sums = [index(track, force=True)[0]['sum'] for track in tracks]
    for left, right in matrix(len(tracks)):
        weight_left, weight_right = normalize(sums[left], sums[right])
        comparisons.append( (left, right, weight_left, weight_right) )

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
                result = metric(x, y)
            results[left, right] = result
        return results

    # Indexed walkers.
    walkers = [walk(track, force_index=True) for track in tracks]

    # Aggregate results.
    totals = defaultdict(lambda: 0)
    counts = defaultdict(lambda: 1)
    for _, _, values in merge(*walkers, merger=merger):
        for comparison, value in values.items():
            if value is not None:
                totals[comparison] += value
                counts[comparison] += 1

    # Create the distance matrix by taking the averages.
    distances = {}
    for comparison in counts:
        distances[comparison] = totals[comparison] / counts[comparison]
    return distances
