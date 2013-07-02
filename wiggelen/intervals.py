"""
Get covered intervals from wiggle tracks and write to BED format.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


import sys


def coverage(walker):
    """
    Get intervals of consecutively defined positions from a walker.

    :arg walker: Tuple of `(region, position, value)` per defined position.
    :type walker: generator(str, int, _)

    :return: Tuples of `(region, begin, end)` per position where `begin` and
        `end` are one-based and inclusive.
    :rtype: generator(str, int, int)

    Example::

        >>> for x in coverage(walk(open('a.wig'))):
        ...     x
        ...
        ('MT', 3, 3)
        ('MT', 5, 20)
        ('MT', 400, 420)
    """
    interval = None

    for region, position, _ in walker:
        if interval is not None:
            if region != interval[0] or position != interval[2] + 1:
                yield interval
                interval = None

        if interval is None:
            interval = region, position, position
        else:
            interval = interval[0], interval[1], position

    # Backlog.
    if interval is not None:
        yield interval


def write(intervals, track=sys.stdout, name=None, description=None):
    """
    Write intervals to a bed track.

    :arg intervals: Tuples of (region, begin, end) per interval.
    :type intervals: generator(str, int, int)
    :arg track: Writable file handle.
    :type track: file
    :arg name: Optional track name (displayed to the left of the track in the
        UCSC Genome Browser).
    :type name: str
    :arg description: Optional track description (displayed as center label in
        the UCSC Genome Browser).
    :type description: str

    Example::

       >>> write(coverage(walk(open('a.wig'))), name='My example')
       track name="My example"
       MT 2 3
       MT 4 20
       MT 399 420
    """
    header = 'track'
    if name is not None:
        header += ' name="%s"' % name
    if description is not None:
        header += ' description="%s"' % description
    header += '\n'
    track.write(header)

    for interval in intervals:
        track.write('%s\t%i\t%i\n' % interval)
