"""
Merge any number of wiggle tracks in various ways.

The algorithm can be parameterized by a merge operation. Four of these
operations are predefined in :attr:`mergers`:

Merger ``sum``: Compute the sum of all values.

Merger ``mean``: Compute the mean of all values (and use 0 for undefined
values).

Merger ``count``: Compute the number of defined values.

Merger ``minus``: Subtract the second value from the first (and use 0 for
undefined values). Only defined on exactly two values.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import division

from .wiggle import zip_


# Compute the sum of all values.
_merger_sum = lambda vs: sum(v for v in vs if v is not None)

# Compute the mean of all values (and use 0 for undefined values).
# Todo: Also provide a mean over only the defined values.
_merger_mean = lambda vs: _merger_sum(vs) / len(vs)

# Compute the number of defined values.
_merger_count = lambda vs: sum(1 for v in vs if v is not None)

# Substract the second from the first (and use 0 for undefined values).
_merger_minus = lambda vs: (vs[0] or 0) - (vs[1] or 0)


#: Predefined mergers. See :mod:`wiggelen.merge` for their definition.
mergers = {'sum':   _merger_sum,
           'mean':  _merger_mean,
           'count': _merger_count,
           'minus': _merger_minus}


def merge(*walkers, **options):
    """
    Merge wiggle tracks.

    This assumes the walkers have their regions in the same order. You can
    force this by using indices. Example::

        >>> from wiggelen import walk
        >>> walkers = [walk(open(track), force_index=True)
        ...            for track in ('a.wig', 'b.wig', 'c.wig')]
        >>> for x in merge(*walkers):
        ...     x
        ...
        ('18', 8, 849.0)
        ('18', 9, 987.0)
        ('MT', 1, 820.0)

    :arg walkers: List of generators yielding tuples of (region, position,
        value) per defined position.
    :type walkers: list(generator(str, int, _))
    :keyword merger: Merge operation (default: sum).
    :type merger: function(list(_) -> _)

    :return: Tuples of (region, position, merged value) per defined position
        in `walkers`.
    :rtype: generator(str, int, _)
    """
    # Todo: Would it be better to also pass region/position to the merger?
    merger = options.get('merger', mergers['sum'])

    for region, position, values in zip_(*walkers):
        yield region, position, merger(values)
