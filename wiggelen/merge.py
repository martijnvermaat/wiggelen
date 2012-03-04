"""
Merge any number of wiggle tracks in various ways.

Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import division

import sys
import argparse

from . import walk_together


# Compute the sum of all values.
_merger_sum = lambda vs: sum(float(v) for v in vs if v is not None)

# Compute the mean of all values (and use 0 for undefined values).
# Todo: Also provide a mean over only the defined values.
_merger_mean = lambda vs: _merger_sum(vs) / len(vs)

# Compute the number of defined values.
_merger_count = lambda vs: sum(1 for v in vs if v is not None)


# Predefined mergers.
mergers = {'sum':   _merger_sum,
           'mean':  _merger_mean,
           'count': _merger_count}


def merge(*walkers, **options):
    """
    Merge wiggle tracks.

    Todo: Is there a better name in (combine, fold, reduce)?
    Todo: Would it be better to also pass region/position to the merger?
    """
    merger = options.get('merger', mergers['sum'])

    for region, position, values in walk_together(*walkers):
        yield region, position, merger(values)
