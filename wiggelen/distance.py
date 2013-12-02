"""
Calculate the distance between two wiggle tracks using a metric designed for
multisets.

This module can be used to assess similarity of next generation sequencing
datasets. A multiset distance measure is used for pairwise comparison of
genomic information as provided by wiggle tracks.

The algorithm can be parameterized by a pairwise distance metric. Four of
these metrics are predefined in :attr:`metrics`:

Metric ``a``: :math:`\\frac{|x - y|}{(x + 1) (y + 1)}`

Metric ``b``: :math:`\\frac{|x - y|}{x + y + 1}`

Metric ``c``: :math:`\\frac{\\text{max}(x, y) \, |x - y|}{(x^2 + 1) (y^2 + 1)}`

Metric ``d``: :math:`\\frac{|x - y|}{\\text{max}(x, y) + 1}`

.. note:: These metrics are ill-defined on the interval (0, 1) so we scale all
       values if necessary.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>
.. moduleauthor:: Jeroen Laros <j.f.j.laros@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import division

from collections import defaultdict
import itertools
import sys

from .wiggle import walk
from .index import Field, index
from .merge import merge


# Todo: Give these metrics a name.
_metric_a = lambda x, y : abs(x - y) / ((x + 1) * (y + 1))
_metric_b = lambda x, y : abs(x - y) / (x + y + 1)
_metric_c = lambda x, y : (max(x, y) * abs(x - y)) \
                          / (((x * x) + 1) * ((y * y) + 1))
_metric_d = lambda x, y : abs(x - y) / (max(x, y) + 1)


#: Predefined pairwise distance metrics. See :mod:`wiggelen.distance` for
#: their definition.
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

    With the default `False` value for `reflexive` and `symmetric`, include
    only the coordinates below the diagonal.

    :arg size: Width and height of the matrix.
    :type size: int
    :arg reflexive: Include coordinates on (`x`, `x`) diagonal.
    :type reflexive: bool
    :arg symmetric: Include coordinates (`x`, `y`) above the diagonal (where
        `x < y`).
    :type symmetric: bool

    :return: All coordinates in the matrix as tuples.
    :rtype: list(int, int)

    Examples::

        >>> matrix(5)
        [(1, 0),
         (2, 0), (2, 1),
         (3, 0), (3, 1), (3, 2),
         (4, 0), (4, 1), (4, 2), (4, 3)]
        >>> matrix(5, reflexive=True)
        [(0, 0),
         (1, 0), (1, 1),
         (2, 0), (2, 1), (2, 2),
         (3, 0), (3, 1), (3, 2), (3, 3),
         (4, 0), (4, 1), (4, 2), (4, 3), (4, 4)]
        >>> matrix(5, symmetric=True)
        [        (0, 1), (0, 2), (0, 3), (0, 4),
         (1, 0),         (1, 2), (1, 3), (1, 4),
         (2, 0), (2, 1),         (2, 3), (2, 4),
         (3, 0), (3, 1), (3, 2),         (3, 4),
         (4, 0), (4, 1), (4, 2), (4, 3)        ]
    """
    return [(i, j) for i in range(0, size) for j in
            itertools.chain(range(i + 1 if reflexive else i),
                            range(i + 1, size if symmetric else 0))]


def distance(*tracks, **options):
    """
    Calculate the pairwise distances between wiggle tracks.

    :arg tracks: List of wiggle tracks.
    :type walkers: list(file)
    :arg metric: Pairwise distance metric (default: a).
    :type merger: function(float, float -> float)
    :arg threshold: Threshold for noise filter (default: no noise filter)
    :type threshold: float

    :return: Pairwise distances between `tracks` as a mapping from
        coordinates in the distance matrix to their values.
    :rtype: dict((int, int), float)

    .. todo:: Check where this goes wrong if we cannot .seek() the tracks.
    .. todo:: Calculate weights per region instead of over the entire track.
    """
    metric = options.get('metric', metrics['a'])
    threshold = options.get('threshold')

    # We construct a list of comparisons for the merger, where each comparison
    # is a tuple of (left, right, weight_left, weight_right).
    comparisons = []

    if threshold:
        field_suffix = '-threshold-%s' % str(threshold)
        noise_filter = lambda value: max(value - threshold, 0)
        def sum_func(acc, value, span):
            return acc + noise_filter(value) * span
        def min_func(acc, value, span):
            filtered = noise_filter(value)
            if filtered > 0:
                return min(acc, noise_filter(value))
            return acc
        fields = [Field('sum' + field_suffix, float, 0, sum_func),
                  Field('posmin' + field_suffix, float, sys.float_info.max,
                        min_func)]
    else:
        field_suffix = ''
        noise_filter = lambda value: value
        fields = []

    summaries = [index(track, force=True, fields=fields)[0]['_all']
                 for track in tracks]

    # Our metrics are undifined on the (0, 1) interval, so if the positive
    # minimum over all tracks is < 1 we upscale everything.
    min_value = min(summary['posmin' + field_suffix] for summary in summaries)
    scale = 1 / min_value if 0 < min_value < 1 else 1

    # Based on the sums of all values in each track we define weights.
    sums = [summary['sum' + field_suffix] for summary in summaries]
    for left, right in matrix(len(tracks)):
        weight_right, weight_left = normalize(sums[left], sums[right])
        comparisons.append( (left, right, weight_left, weight_right) )

    # The merger returns a matrix of values, one for each pairwise comparison.
    def merger(values):
        results = {}
        for left, right, weight_left, weight_right in comparisons:
            value_left, value_right = values[left], values[right]
            if value_left is None and value_right is None:
                result = None
            else:
                x = (weight_left * scale * noise_filter(value_left)
                     if value_left else 0)
                y = (weight_right * scale * noise_filter(value_right)
                     if value_right else 0)
                result = metric(x, y) if x or y else None
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
