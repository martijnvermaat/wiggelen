"""
Read and write wiggle tracks.

Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


import sys
import itertools


def walk(track=sys.stdin):
    """
    Walk over the track and for each position yield the position and value.

    @arg track: Wiggle track.
    @type track: file

    @return: Pairs of (position, value) per defined position.
    @rtype: generator(int, str)

    Todo: Handle chromosomes.
    Todo: Handle variable and fixed step.
    """
    for line in track:
        fields = line.split()
        try:
            position = int(fields[0])
            # Todo: Make value of int type if it is an integer
            value = float(fields[1])
        except ValueError:
            continue
        yield position, value


def walk_together(*walkers):
    """
    Walk over all tracks simultaneously and for each position yield the
    position and a list of values for each track, or None in case the track
    has no value on the position.

    @arg walkers: List of generators yielding pairs of (position, value) per
        defined position.
    @rtype: list(generator(int, str))

    @return: Pairs of (position, values) per defined position.
    @rtype: generator(int, list(str))

    Todo: Handle chromosomes and their order.
    """
    items = []
    for walker in walkers:
        try:
            items.append(next(walker))
        except StopIteration:
            items.append(None)

    while True:
        if not any(items):
            break

        position = min(item[0] for item in items if item is not None)

        values = [item[1] if item is not None and item[0] == position else None
                  for item in items]
        yield position, values

        for i, item in enumerate(items):

            if item is not None and item[0] == position:
                try:
                    items[i] = next(walkers[i])
                except StopIteration:
                    items[i] = None


# Filter items from a walker
filter_ = itertools.ifilter


def write(walker, track=sys.stdout):
    """
    Write items from a walker to a wiggle track.

    Todo: Handle chromosomes.
    Todo: Options for variable or fixed step, window size, etc.
    """
    track.write('track type=wiggle_0 name="" description=""\n')
    for item in walker:
        # Todo: See if float/int types both are represented correctly
        track.write('%d %s\n' % item)
