"""
Read and write wiggle tracks.

The `wiggle (WIG) format <https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html>`_
is for storing dense, continuous genomic data such as GC percent, probability
scores, read depth, and transcriptome data.

.. todo: Note in the documentation that walker values can be of any type, but
    that valid wiggle tracks only have int or float values.
.. todo: Integrate some of our existing scripts: ``ngs-misc/sageWiggle``,
    ``gapss3/pileup2wig``, ``gapss3/mpileup2wig``, ``ngs-data/wiggle2region``.
.. todo: Note that positioning is one-based.
.. todo: Note in manual that itertools.ifilter and itertools.imap are useful.

.. Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
.. Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
.. Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
"""


import sys

from .index import index, write_index


def walk(track=sys.stdin, force_index=False):
    """
    Walk over the track and yield (region, position, value) tuples.

    The values are always of type `int` or `float`.

    :arg track: Wiggle track.
    :type track: file
    :arg force_index: Force creating an index if it does not yet exist.
    :type force_index: bool

    :return: Tuples of (region, position, value) per defined position.
    :rtype: generator(str, int, _)

    Example::

        >>> for x in walk():
        ...     x
        ('chr18', 34344, 629.0)
        ('chr18', 34345, 649.0)
        ('chr18', 34446, 657.0)
        ('chrM',  308,   520.0)
        ('chrM',  309,   519.0)

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
        _, mapping, _ = idx
        # Todo: Sort in a way that is compatible with existing wiggle tracks.
        #     Inspiration could be sorted BAM files. GATK requires these to be
        #     sorted according to the order in the reference file.
        regions = sorted(mapping)

    for expected_region in regions:

        if expected_region is not None:
            track.seek(mapping[expected_region])

        for line in track:

            if line.startswith('browser') or line.startswith('track') \
                   or line.startswith('#') or line.isspace():
                # As far as I can see empty lines and comments are not allowed
                # by the spec, but I guess they exist in the real world.
                pass

            elif line.startswith('variableStep'):
                try:
                    fields = dict(map(lambda field: field.split('='),
                                      line[len('variableStep'):].split()))
                    region = fields['chrom']
                    span = fields.get('span', 1)
                    format_ = 'variable'
                except (ValueError, KeyError):
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
                except (ValueError, KeyError):
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

    # Todo: Automatically build the index (if we read the whole file) and
    # write it to a file?


def zip_(*walkers):
    """
    Walk over all tracks simultaneously and for each position yield the
    region, position and a list of values for each track, or ``None`` in case
    the track has no value on the position.

    .. note:: This assumes the order of regions is compatible over all
        walkers. If you are unsure if this is the case for your input wiggle
        tracks, use the :func:`walk` function with the ``force_index`` keyword
        argument.

    :arg walkers: List of generators yielding tuples of (region, position,
        value) per defined position.
    :type walkers: list(generator(str, int, _))

    :return: Tuples of (region, position, values) per defined position.
    :rtype: generator(str, int, list(_))

    Example::

        >>> for x in zip_(walk(open('a.wig')), walk(open('b.wig'))):
        ...     x
        ('18', 7, [29.0, None])
        ('18', 8, [49.0, None])
        ('18', 9, [None, 87.0])
        ('MT', 1, [20.0, None])
        ('MT', 2, [36.0, 92.0])
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


def fill(walker, regions=None):
    """
    Fill in undefined positions with None values.

    :arg walker: Tuple of (region, position, value) per defined position.
    :type walker: generator(str, int, _)
    :arg regions: Dictionary with regions as keys and (start, stop) tuples as
        values. If not ``None``, fill positions from start to stop (both
        including) in these regions. If ``None``, fill positions in all
        regions between their first and last defined positions.
    :type regions: dict(str, (int, int))

    :return: Tuples of (region, position, value) per position where value is
        ``None`` if it was not defined in the original walker.
    :rtype: generator(str, int, _)

    Example::

        >>> for x in walk(open('a.wig')):
        ...     x
        ('MT', 3, 29.0)
        ('MT', 5, 49.0)
        ('MT', 8, 87.0)
        ('MT', 9, 20.0)
        >>> for x in fill(walk(open('a.wig')):
        ...     x
        ('MT', 3, 29.0)
        ('MT', 4, None)
        ('MT', 5, 49.0)
        ('MT', 6, None)
        ('MT', 7, None)
        ('MT', 8, 87.0)
        ('MT', 9, 20.0)
    """
    previous_region = previous_position = None

    for region, position, value in walker:
        if region != previous_region:
            # Backlog.
            if regions is not None:
                try:
                    start, stop = regions[previous_region]
                    for p in xrange(max(previous_position + 1, start), stop + 1):
                        yield previous_region, p, None
                except KeyError:
                    pass
            previous_region = region
            previous_position = None

        if regions is None:
            # No explicitely specified regions, fill everything.
            if previous_position is not None:
                for p in xrange(previous_position + 1, position):
                    yield region, p, None
        else:
            try:
                # Specified where we must fill.
                start, stop = regions[region]
                if previous_position is None:
                    for p in xrange(start, min(position, stop + 1)):
                        yield region, p, None
                else:
                    for p in xrange(max(previous_position + 1, start),
                                    min(position, stop + 1)):
                        yield region, p, None
            except KeyError:
                # Region is not in explicitely specified regions, don't
                # fill it.
                pass

        previous_position = position
        yield region, position, value

    # Backlog.
    if regions is not None:
        try:
            start, stop = regions[previous_region]
            for p in xrange(max(previous_position + 1, start), stop + 1):
                yield previous_region, p, None
        except KeyError:
            pass


def write(walker, track=sys.stdout, serializer=str):
    """
    Write items from a walker to a wiggle track.

    :arg walker: Tuples of (region, position, value) per defined position.
    :type walker: generator(str, int, _)
    :arg track: Writable file handle.
    :type track: file
    :arg serializer: Function making strings from values.
    :type serializer: function(_ -> str)

    Example::

       >>> write(walk(open('a.wig')))
       track type=wiggle_0 name="" description=""
       variableStep chrom=1
       1 520.0
       4 536.0
       8 553.0
       variableStep chrom=MT
       1 568.0
       2 598.0
       6 616.0

    .. todo:: Options for variable or fixed step, window size, etc.
    """
    size = 0

    header = 'track type=wiggle_0 name="" description=""\n'
    track.write(header)
    size += len(header)

    summary = {'sum': 0, 'count': 0}
    mapping = {}
    current_region = None

    for region, position, value in walker:
        if region != current_region:
            mapping[region] = size
            chrom = 'variableStep chrom=%s\n' % region
            track.write(chrom)
            size += len(chrom)
            current_region = region
        step = '%d %s\n' % (position, serializer(value))
        track.write(step)
        size += len(step)
        summary['sum'] += value
        summary['count'] += 1

    write_index(summary, mapping, track)
