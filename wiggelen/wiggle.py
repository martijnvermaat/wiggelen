"""
Read and write wiggle tracks.

.. todo:: Now that we also include the end positions in the index, it is
    possible to do random jumps inside a region with some educated guessing.
    Perfect hits would not be possible, since the length of the lines is
    variable.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


import sys

from .parse import LineType, create_state, parse
from .index import ReadError, index, write_index


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
        ...
        ('chr18', 34344, 629.0)
        ('chr18', 34345, 649.0)
        ('chr18', 34446, 657.0)
        ('chrM',  308,   520.0)
        ('chrM',  309,   519.0)

    .. todo:: Optionally give a list of regions to walk, in that order.
    """
    # Todo: Do something with browser and track lines.
    # Todo: Better exceptions.
    # Todo: Prettify this code.
    # Todo: Detect if index does not agree with track.
    region = None

    idx, _ = index(track, force=force_index)

    if idx is None:
        regions = [None]
    else:
        # Todo: Sort in a way that is compatible with existing wiggle tracks.
        #     Inspiration could be sorted BAM files. GATK requires these to be
        #     sorted according to the order in the reference file.
        regions = sorted(r for r in idx if r != '_all')

    for expected_region in regions:

        if expected_region is not None:
            track.seek(idx[expected_region]['start'])

        state = create_state()

        for line in track:
            line_type, data = parse(line, state)
            if line_type == LineType.REGION:
                region = data
                if expected_region is not None and region != expected_region:
                    break
            elif line_type == LineType.DATA:
                # Optimization: A `while` loop is faster than `for` and
                # `range`.
                i = 0
                while i < data.span:
                    yield region, data.position + i, data.value
                    i += 1

        # Todo: If there is no index yet, but we read the whole file, write
        #     the index we just built. But if we go this route, I think we
        #     should somehow merge wiggle.walk with index.index since they
        #     will be doing virtually the same thing.
        #
        #     if expected_region is None:
        #        write_index(idx, track)


def zip_(*walkers):
    """
    Walk over all tracks simultaneously and for each position yield the
    region, position and a list of values for each track, or `None` in case
    the track has no value on the position.

    .. note:: This assumes the order of regions is compatible over all
        walkers. If you are unsure if this is the case for your input wiggle
        tracks, use the :func:`walk` function with the `force_index` keyword
        argument.

    :arg walkers: List of generators yielding tuples of (region, position,
        value) per defined position.
    :type walkers: list(generator(str, int, _))

    :return: Tuples of (region, position, values) per defined position.
    :rtype: generator(str, int, list(_))

    Example::

        >>> for x in zip_(walk(open('a.wig')), walk(open('b.wig'))):
        ...     x
        ...
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


def fill(walker, regions=None, filler=None, only_edges=False):
    """
    Fill in undefined positions with `filler` (or `None`).

    :arg walker: Tuple of (region, position, value) per defined position.
    :type walker: generator(str, int, _)
    :arg regions: Dictionary with regions as keys and (start, stop) tuples as
        values. If not `None`, fill positions from start to stop (both
        including) in these regions. If `None`, fill positions in all
        regions between their first and last defined positions.
    :type regions: dict(str, (int, int))
    :arg filler: Value to use for filling undefined positions.
    :type filler: _
    :arg only_edges: Only fill the first and last of continuously undefined
        positions.
    :type only_edges: bool

    :return: Tuples of (region, position, value) per position where value is
        `filler` if it was not defined in the original walker.
    :rtype: generator(str, int, _)

    Example::

        >>> for x in walk(open('a.wig')):
        ...     x
        ...
        ('MT', 3, 29.0)
        ('MT', 5, 49.0)
        ('MT', 8, 87.0)
        ('MT', 9, 20.0)
        >>> for x in fill(walk(open('a.wig')):
        ...     x
        ...
        ('MT', 3, 29.0)
        ('MT', 4, None)
        ('MT', 5, 49.0)
        ('MT', 6, None)
        ('MT', 7, None)
        ('MT', 8, 87.0)
        ('MT', 9, 20.0)

    The `only_edges` argument might seem a bit out of place here, but can be
    useful in combination with `filler=0` when creating a line plot. Without
    any filling, non-zero lines may be plotted where there is actually no
    data.

    .. note:: This might be a tiny bit memory-hungry on Python 2.x if there
        are *very* large gaps to fill since we use the range function to
        generate the positions. I don't think it's worth it to add version
        specific code paths for this.
    """
    previous_region = previous_position = None

    for region, position, value in walker:
        if region != previous_region:
            # Backlog.
            if regions is not None:
                try:
                    start, stop = regions[previous_region]
                    step = max((stop - max(previous_position + 1, start))
                               * only_edges, 1)
                    for p in range(max(previous_position + 1, start),
                                   stop + 1, step):
                        yield previous_region, p, filler
                except KeyError:
                    pass
            previous_region = region
            previous_position = None

        if regions is None:
            # No explicitely specified regions, fill everything.
            if previous_position is not None:
                step = max((position - previous_position - 2) * only_edges,
                           1)
                for p in range(previous_position + 1, position, step):
                    yield region, p, filler
        else:
            try:
                # Specified where we must fill.
                start, stop = regions[region]
                if previous_position is None:
                    step = max((min(position, stop + 1) - start - 1)
                               * only_edges, 1)
                    for p in range(start, min(position, stop + 1), step):
                        yield region, p, filler
                else:
                    step = max((min(position, stop + 1) -
                                max(previous_position + 1, start) - 1)
                               * only_edges, 1)
                    for p in range(max(previous_position + 1, start),
                                   min(position, stop + 1), step):
                        yield region, p, filler
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
            step = max((stop - max(previous_position + 1, start))
                       * only_edges, 1)
            for p in range(max(previous_position + 1, start), stop + 1, step):
                yield previous_region, p, filler
        except KeyError:
            pass


def write(walker, track=sys.stdout, serializer=str, name=None,
          description=None):
    """
    Write items from a walker to a wiggle track.

    :arg walker: Tuples of (region, position, value) per defined position.
    :type walker: generator(str, int, _)
    :arg track: Writable file handle.
    :type track: file
    :arg serializer: Function making strings from values.
    :type serializer: function(_ -> str)
    :arg name: Optional track name (displayed to the left of the track in the
        UCSC Genome Browser).
    :type name: str
    :arg description: Optional track description (displayed as center label in
        the UCSC Genome Browser).
    :type description: str

    .. note:: Values of `None` are discarded.

    Example::

       >>> write(walk(open('a.wig')), name='My example')
       track type=wiggle_0 name="My example"
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

    header = 'track type=wiggle_0'
    if name is not None:
        header += ' name="%s"' % name
    if description is not None:
        header += ' description="%s"' % description
    header += '\n'
    track.write(header)
    size += len(header)

    idx = {}
    current_region = None

    for region, position, value in walker:
        if value is None:
            continue
        if region != current_region:
            line = 'variableStep chrom=%s\n' % region
            track.write(line)
            idx[region] = {
                'region': region,
                'start':  size,
                'stop':   size + len(line),
                'sum':    0,
                'min':    sys.float_info.max,
                'posmin': sys.float_info.max,
                'max':    0,
                'count':  0}
            size += len(line)
            current_region = region
        line = '%d %s\n' % (position, serializer(value))
        track.write(line)
        idx[region]['stop'] = size + len(line)
        idx[region]['sum'] += value
        idx[region]['min'] = min(value, idx[region]['min'])
        if value > 0:
            idx[region]['posmin'] = min(value, idx[region]['posmin'])
        idx[region]['max'] = max(value, idx[region]['max'])
        idx[region]['count'] += 1
        size += len(line)

    idx['_all'] = {
        'region': '_all',
        'start':  0,
        'stop':   size,
        'sum':    sum(r['sum'] for r in idx.values()),
        'min':    min(r['min'] for r in idx.values()),
        'posmin': min(r['posmin'] for r in idx.values()),
        'max':    max(r['max'] for r in idx.values()),
        'count':  sum(r['count'] for r in idx.values())}

    write_index(idx, track)
