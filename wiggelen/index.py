"""
Index regions/chromosomes in wiggle tracks for random access.

Indexing a wiggle track results in a mapping of regions to summaries. The
summaries are dictionaries including start and stop positions and some
statistical metrics. A summary for the entire wiggle track is included as the
special region ``_all``.

This data can be written to a file next to the wiggle track file (in case this
is a regular file). Example of the serialization we use::

    region=_all,start=0,stop=12453,sum=4544353,count=63343
    region=1,start=47,stop=3433,sum=4353,count=643
    region=X,start=3433,stop=8743,sum=454,count=343
    region=Y,start=8743,stop=10362,sum=7353,count=343
    region=MT,start=10362,stop=12453,sum=353,count=143

Note that we do not impose a certain order on the lines in the index nor on
the fields on a line.

.. todo:: Add some other metrics to the index (standard deviation, min, max).

.. todo:: Mechanism to include custom metrics in the index?

.. todo:: Cache the index objects somehow during the process. Unfortunately,
    we cannot attach it to the track file handler, as it does not accept
    additional attributes. I guess we can do something based on the hash of
    the file handler?

.. Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
.. Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
.. Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
"""


import sys
from collections import defaultdict

from .parse import LineType, create_state, parse


#: Whether or not indices are written to a file.
WRITE_INDEX = True

#: Suffix used for index files.
INDEX_SUFFIX = '.idx'


def _cast(field, value):
    casters = defaultdict(lambda: str,
                          start=int, stop=int, sum=float, count=int)
    return casters[field](value)


def _index_filename(track=sys.stdin):
    """
    Try to create a filename for the index file.
    """
    filename = getattr(track, 'name', None)
    if filename is not None and not filename.startswith('<'):
        return filename + INDEX_SUFFIX


def write_index(idx, track=sys.stdout):
    """
    Try to write the index to a file and return its filename.

    :arg idx: Wiggle track index.
    :type idx: dict(str, dict(str, _))
    :arg track: Wiggle track the index belongs to.
    :type track: file

    :return: Filename for the written index, or ``None`` if the index could
        not be written.
    :rtype: str
    """
    if not WRITE_INDEX:
        return

    filename = _index_filename(track)

    if filename is None:
        return

    try:
        with open(filename, 'w') as f:
            f.write('\n'.join(','.join('%s=%s' % d for d in s.items())
                              for s in idx.values()) + '\n')
        return filename
    except IOError:
        pass


def read_index(track=sys.stdin):
    """
    Try to read the index from a file.

    :arg track: Wiggle track the index belongs to.
    :type track: file

    :return: Wiggle track index, or ``None`` if the index could not be read.
    :rtype: dict(str, dict(str, _))

    .. todo:: Only accept if index is newer than wiggle track?
    """
    filename = _index_filename(track)

    if filename is None:
        return

    try:
        with open(filename) as f:
            idx = {}
            for line in f:
                summary = {}
                for d in line.rstrip().split(','):
                    k, v = d.split('=')
                    summary[k] = _cast(k, v)
                idx[summary['region']] = summary
            return idx
    except IOError:
        pass


def index(track=sys.stdin, force=False):
    """
    Return index of region positions in track.

    :arg track: Wiggle track.
    :type track: file
    :arg force: Force creating an index if it does not yet exist.
    :type force: bool

    :return: Wiggle track index and index filename.
    :rtype: dict(str, dict(str, _)), str

    .. todo:: Handle non-writable index, corrupt index, etc.
    """
    idx = read_index(track)

    if idx is not None or not force:
        # Todo: Shouldn't we return a filename here?
        return idx, None

    region = None
    idx = {}
    state = create_state()

    while True:
        line = track.readline()
        if not line:
            break
        line_type, data = parse(line, state)

        if line_type == LineType.REGION:
            region = data
            idx[region] = {
                'region': region,
                'start':  track.tell() - len(line),
                'stop':   track.tell(),
                'sum':    0,
                'count':  0}
        elif line_type == LineType.DATA:
            idx[region]['stop'] = track.tell()
            idx[region]['sum'] += data.value * data.span
            idx[region]['count'] += data.span

    # Todo: Technically, start of _all is not 0 since we always have a track
    # line in the file.
    idx['_all'] = {
        'region': '_all',
        'start':  0,
        'stop':   track.tell(),
        'sum':    sum(r['sum'] for r in idx.values()),
        'count':  sum(r['count'] for r in idx.values())}

    return idx, write_index(idx, track)
