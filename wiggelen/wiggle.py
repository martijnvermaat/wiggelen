"""
Read and write wiggle tracks.

Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


import sys
import itertools

from .index import index


def walk(track=sys.stdin, regions=None, index=None):
    """
    Walk over the track and yield (region, position, value) triples.

    @arg track: Wiggle track.
    @type track: file
    @arg regions: Optional list of regions to walk, in that order.
    @type regions: list(str)
    @arg index: Optional index to use for finding regions.
    @type index: dict(str, int)

    @return: Triples of (region, position, value) per defined position.
    @rtype: generator(str, int, str)

    Todo: Do something with browser and track lines.
    Todo: Better exceptions.
    Todo: Prettify the parsing code.
    Todo: By default walk the regions from the index in alphabetical order,
        even if no explicit list of regions is given. This would remove the
        need for ordered_regions with merge.
    """
    format = region = start = step = span = None

    #Todo: Work with the supplied regions and index.
    #for region in regions:
    #    if index:
    #        track.seek(index[region])
    #    else:
    #        # check if we get expected region (possibly skipping empty ones)

    for line in track:
        if line.startswith('browser') or line.startswith('track'):
            pass

        elif line.startswith('variableStep'):
            try:
                fields = dict(map(lambda field: field.split('='),
                                  line[len('variableStep'):].split()))
                region = fields['chrom']
                span = fields.get('span', 1)
                format = 'variable'
            except ValueError, KeyError:
                raise Exception('Could not parse line: %s' % line)

        elif line.startswith('fixedStep'):
            try:
                fields = dict(map(lambda field: field.split('='),
                                  line[len('fixedStep'):].split()))
                region = fields['chrom']
                start = int(fields['start'])
                step = int(fields['step'])
                span = int(fields.get('span', 1))
                format = 'fixed'
            except ValueError, KeyError:
                raise Exception('Could not parse line: %s' % line)

        elif format == 'variable':
            try:
                position, value = line.split()
                position = int(position)
                value = float(value) if '.' in value else int(value)
            except ValueError:
                raise Exception('Could not parse line: %s' % line)
            for i in range(span):
                yield region, position + i, value

        elif format == 'fixed':
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
    region, position and a list of values for each track, or None in case the
    track has no value on the position.

    @arg walkers: List of generators yielding triples of (region, position,
        value) per defined position.
    @rtype: list(generator(str, int, str))

    @return: Triples of (region, position, values) per defined position.
    @rtype: generator(str, int, list(str))
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

        region, position = min(item[0:2] for item in items if item is not None)

        values = [item[1] if item is not None and item[0:2] == [region, position]
                  else None for item in items]
        yield region, position, values

        for i, item in enumerate(items):

            if item is not None and item[0:2] == [region, position]:
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
    for region, position, value in walker:
        # Todo: See if float/int types both are represented correctly
        track.write('%d %s\n' % (position, value))
