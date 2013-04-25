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

Additional custom fields can be added to the index by providing custom field
definitions. Such a definition is a 4-tuple with the following fields:

* The name of the field.
* A function casting a field value from `string`.
* Initial value.
* Reducer function taking as inputs the accumulated field value, the current
  value, and the current span, returning a new accumulated field value.

As an example, the standard `sum` field could be defined as the
following tuple::

    'sum', float, 0, lambda acc, value, span: acc + value * span

In practice, choose unique names for custom fields.

.. todo:: Add some other metrics to the index (standard deviation, min, max).

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


import sys
from collections import defaultdict

from .parse import LineType, create_state, parse


#: Whether or not indices are written to a file.
WRITE_INDEX = True

#: Suffix used for index files.
INDEX_SUFFIX = '.idx'

#: Whether or not indices are cached in memory during execution.
CACHE_INDEX = True


# Cache store of indices, indexed by index filename.
_cache = {}


class ReadError(Exception):
    """
    Raised if a wiggle track does not provide random access. Reading with
    random access is needed for using and creating an index.
    """
    pass


# Cast index values to their respective types.
def _cast(field, value, customs=None):
    customs = customs or []
    casters = defaultdict(lambda: str,
                          start=int, stop=int, sum=float, count=int)
    casters.update(dict((field, caster) for field, caster, _, _ in customs))
    return casters[field](value)


# Try to create a filename for the index file.
def _index_filename(track=sys.stdin):
    filename = getattr(track, 'name', None)
    if filename is not None and not filename.startswith('<'):
        # Note that this does not play nice if we have changed our current
        # working dir and the file was opened using a relative path.
        return filename + INDEX_SUFFIX


def clear_cache():
    """
    Clear the in-memory cache of index objects.
    """
    _cache.clear()


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
    filename = _index_filename(track)

    if filename is None:
        return

    if CACHE_INDEX:
        _cache[filename] = idx

    if not WRITE_INDEX:
        return

    try:
        with open(filename, 'w') as f:
            f.write('\n'.join(','.join('%s=%s' % d for d in s.items())
                              for s in idx.values()) + '\n')
        return filename
    except IOError:
        pass


def read_index(track=sys.stdin, customs=None):
    """
    Try to read the index from a file.

    :arg track: Wiggle track the index belongs to.
    :type track: file
    :arg customs: List of custom index field definitions.
    :type customs: list

    :return: Wiggle track index, or ``None`` if the index could not be read.
    :rtype: dict(str, dict(str, _))

    .. todo:: Only accept if index is newer than wiggle track?
    """
    customs = customs or []

    filename = _index_filename(track)

    if filename is None:
        return

    if CACHE_INDEX and filename in _cache:
        return _cache[filename]

    try:
        idx = {}
        with open(filename) as f:
            for line in f:
                summary = {}
                for d in line.rstrip().split(','):
                    k, v = d.split('=')
                    summary[k] = _cast(k, v, customs=customs)
                idx[summary['region']] = summary
        # Todo: Here we check if all the required custom fields are in the
        #     index for `_all`, but we should really do a better check if all
        #     required data is there.
        if all(field in idx['_all'] for field, _, _, _ in customs):
            return idx
    except IOError:
        pass


def index(track=sys.stdin, force=False, customs=None):
    """
    Return index of region positions in track.

    :arg track: Wiggle track.
    :type track: file
    :arg force: Force creating an index if it does not yet exist.
    :type force: bool
    :arg customs: List of custom index field definitions.
    :type customs: list

    :return: Wiggle track index and index filename.
    :rtype: dict(str, dict(str, _)), str

    .. todo:: It is not possible to force the index to be rewritten if it
        already exists.
    .. todo:: Handle non-writable index, corrupt index, etc.
    """
    customs = customs or []

    idx = read_index(track, customs=customs)

    if idx is not None or not force:
        return idx, _index_filename(track)

    try:
        track.tell()
    except IOError:
        raise ReadError('Could not index track (needs random access)')

    region = None
    idx = {'_all': {'region': '_all',
                    'start':  0,
                    'stop':   track.tell(),
                    'sum':    0,
                    'count':  0}}
    idx['_all'].update(dict((field, init) for field, _, init, _ in customs))

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
            idx[region].update(dict((field, init)
                                    for field, _, init, _ in customs))
        elif line_type == LineType.DATA:
            for r in region, '_all':
                idx[r]['stop'] = track.tell()
                idx[r]['sum'] += data.value * data.span
                idx[r]['count'] += data.span
                idx[r].update(dict((field, reducer(idx[r][field],
                                                   data.value,
                                                   data.span))
                               for field, _, _, reducer in customs))

    return idx, write_index(idx, track)
