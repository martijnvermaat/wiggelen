"""
Read and write wiggle tracks.

.. Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
.. Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
.. Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
"""


import sys
import itertools

from .index import index, write_index


def walk(track=sys.stdin, force_index=False):
    """
    Walk over the track and yield (region, position, value) triples.

    :arg track: Wiggle track.
    :type track: file
    :arg force_index: Force creating an index if it does not yet exist.
    :type force_index: bool

    :return: Triples of (region, position, value) per defined position.
    :rtype: generator(str, int, float)

    .. todo:: Optionally give a list of regions to walk, in that order.
    .. todo:: Do something with browser and track lines.
    .. todo:: Better exceptions.
    .. todo:: Prettify this code.
    .. todo:: Detect if index does not agree with track.
    """
    format_ = region = start = step = span = None

    idx = index(track, force=force_index)

    if idx is None:
        regions = [None]
    else:
        # Todo: Sort in a way that is compatible with existing wiggle tracks.
        #     Inspiration could be sorted BAM files. GATK requires these to be
        #     sorted according to the order in the reference file.
        regions = sorted(idx)

    for expected_region in regions:

        if expected_region is not None:
            track.seek(idx[expected_region])

        for line in track:

            if line.startswith('browser') or line.startswith('track'):
                pass

            elif line.startswith('variableStep'):
                try:
                    fields = dict(map(lambda field: field.split('='),
                                      line[len('variableStep'):].split()))
                    region = fields['chrom']
                    span = fields.get('span', 1)
                    format_ = 'variable'
                except ValueError, KeyError:
                    raise Exception('Could not parse line: %s' % line)
                if expected_region is not None and region != expected_region:
                    break

            elif line.startswith('fixedStep'):
                try:
                    fields = dict(map(lambda field: field.split('='),
                                      line[len('fixedStep'):].split()))
                    region = fields['chrom']
                    start = int(fields['start'])
                    step = int(fields['step'])
                    span = int(fields.get('span', 1))
                    format_ = 'fixed'
                except ValueError, KeyError:
                    raise Exception('Could not parse line: %s' % line)
                if expected_region is not None and region != expected_region:
                    break

            elif format_ == 'variable':
                try:
                    position, value = line.split()
                    position = int(position)
                    value = float(value) if '.' in value else int(value)
                except ValueError:
                    raise Exception('Could not parse line: %s' % line)
                for i in range(span):
                    yield region, position + i, value

            elif format_ == 'fixed':
                try:
                    value = float(line) if '.' in line else int(line)
                except ValueError:
                    raise Exception('Could not parse line: %s' % line)
                for i in range(span):
                    yield region, start + i, value
                start += step

            else:
                raise Exception('Could not parse line: %s' % line)


def walk_together(*walkers):
    """
    Walk over all tracks simultaneously and for each position yield the
    region, position and a list of values for each track, or ``None`` in case
    the track has no value on the position.

    .. note:: This assumes the order of regions is compatible over all
        walkers. If you are unsure if this is the case for your input wiggle
        tracks, use the :func:`walk` function with the ``force_index`` keyword
        argument.

    :arg walkers: List of generators yielding triples of (region, position,
        value) per defined position.
    :type walkers: list(generator(str, int, float))

    :return: Triples of (region, position, values) per defined position.
    :rtype: generator(str, int, list(float))
    """
    # We work with a list of lookahead items. If a walker has no more items,
    # we use None in the lookahead list.
    items = []
    for walker in walkers:
        try:
            items.append(next(walker))
        except StopIteration:
            items.append(None)

    # Regions seen so far.
    regions = set()
    previous_region = None

    while True:
        # If all lookahead items are None, we are done.
        if not any(items):
            break

        # Get the next position to yield.
        region, position = min(item[0:2] for item in items if item is not None)

        # Check region order compatibility.
        if region != previous_region:
            if region in regions:
                raise Exception('The order of regions is not compatible')
            regions.add(region)
            previous_region = region

        # Yield all values at this position.
        values = [item[2] if item is not None and item[0:2] == (region, position)
                  else None for item in items]
        yield region, position, values

        # Advance the lookahead list where we just yielded a value.
        for i, item in enumerate(items):
            if item is not None and item[0:2] == (region, position):
                try:
                    items[i] = next(walkers[i])
                except StopIteration:
                    items[i] = None


#: Filter items from a walker.
filter_ = itertools.ifilter


def write(walker, track=sys.stdout):
    """
    Write items from a walker to a wiggle track.

    :arg walker: Triples of (region, position, value) per defined position.
    :type walker: generator(str, int, float)

    .. todo:: Options for variable or fixed step, window size, etc.
    """
    size = 0

    header = 'track type=wiggle_0 name="" description=""\n'
    track.write(header)
    size += len(header)

    idx = {}
    current_region = None

    for region, position, value in walker:
        if region != current_region:
            idx[region] = size
            chrom = 'variableStep chrom=%s\n' % region
            track.write(chrom)
            size += len(chrom)
            current_region = region
        # Todo: See if float/int types both are represented correctly.
        step = '%d %s\n' % (position, value)
        track.write(step)
        size += len(step)

    write_index(idx, track)
