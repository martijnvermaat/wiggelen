"""
Index regions/chromosomes in wiggle tracks for random access.

Indexing a wiggle track results in two dictionaries: the first contains some
summary data, the second is a mapping of regions and their positions in the
wiggle track file.

This data can be written to a file next to the wiggle track file (in case this
is a regular file). Example of the serialization we use:

    #sum=4544353,count=63343
    1 47
    X 3433
    Y 8743
    MT 10362

.. todo:: Cache the index objects somehow during the process. Unfortunately,
    we cannot attach it to the track file handler, as it does not accept
    additional attributes.

.. Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
.. Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
.. Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

.. Licensed under the MIT license, see the LICENSE file.
"""


import sys


#: Suffix used for index files.
INDEX_SUFFIX = '.idx'


def _index_filename(track=sys.stdin):
    """
    Try to create a filename for the index file.
    """
    filename = getattr(track, 'name', None)
    if filename is not None and not filename.startswith('<'):
        return filename + INDEX_SUFFIX


def write_index(summary, mapping, track=sys.stdout):
    """
    Try to write the index to a file and return its filename.

    :arg summary: Wiggle track summary data.
    :type summary: dict(str, _)
    :arg mapping: Mapping of regions to wiggle track positions.
    :type mapping: dict(str, int)
    :arg track: Wiggle track the index belongs to.
    :type track: file

    :return: Filename for the written index, or ``None`` if the index could
        not be written.
    :rtype: str
    """
    filename = _index_filename(track)

    if filename is not None:
        try:
            with open(filename, 'w') as f:
                f.write('#' + ','.join('%s=%s' % d for d in summary.items()) + '\n')
                f.write('\n'.join('%s %d' % i for i in mapping.items()) + '\n')
                return filename
        except IOError:
            pass


def read_index(track=sys.stdin):
    """
    Try to read the index from a file.

    :arg track: Wiggle track the index belongs to.
    :type track: file

    :return: Wiggle track summary and mapping, or ``None`` if the index could
        not be read.
    :rtype: tuple(dict(str, _), dict(str, int))

    .. todo:: Only accept if index is newer than wiggle track?
    """
    filename = _index_filename(track)

    if filename is not None:
        try:
            with open(filename) as f:
                summary = dict((k, float(v)) for k, v in
                               (d.split('=') for d in next(f)[1:-1].split(',')))
                mapping = dict((r, int(p)) for r, p in (l.split() for l in f))
                return summary, mapping
        except IOError:
            pass


def index(track=sys.stdin, force=False):
    """
    Return index of region positions in track.

    :arg track: Wiggle track.
    :type track: file
    :arg force: Force creating an index if it does not yet exist.
    :type force: bool

    :return: Wiggle track summary and mapping.
    :rtype: tuple(dict(str, _), dict(str, int))

    .. todo:: Handle non-writable index, corrupt index, etc.
    .. todo:: Also including the end positions would make it possible to do
        random jumps inside a region with some educated guessing. Perfect hits
        would not be possible, since the length of the lines are variable.
    .. todo:: Also return index filename so we can check if it was written.
    .. todo:: Just to calculate the summaries, we repeat ourselves quite a bit
        here compared to the parsing in wiggle.walk, we could refactor this.
    """
    idx = read_index(track)

    if idx is not None:
        summary, mapping = idx
        return summary, mapping, None

    if not force:
        return None

    # Todo: Populate summary.
    summary = {'sum': 0, 'count': 0}
    mapping = {}
    format_ = span = None
    while True:
        line = track.readline()
        if not line:
            break
        if line.startswith('browser') or line.startswith('track'):
            pass
        elif line.startswith('variableStep'):
            region = line.split('chrom=')[1].split()[0]
            mapping[region] = track.tell() - len(line)
            try:
                span = line.split('span=')[1].split()[0]
            except IndexError:
                span = 1
            format_ = 'variable'
        elif line.startswith('fixedStep'):
            region = line.split('chrom=')[1].split()[0]
            mapping[region] = track.tell() - len(line)
            try:
                span = line.split('span=')[1].split()[0]
            except IndexError:
                span = 1
            format_ = 'fixed'
        elif format_ == 'variable':
            try:
                position, value = line.split()
                value = float(value) if '.' in value else int(value)
            except ValueError:
                raise Exception('Could not parse line: %s' % line)
            summary['count'] += 1 * span
            summary['sum'] += value * span
        elif format_ == 'fixed':
            try:
                value = float(line) if '.' in line else int(line)
            except ValueError:
                raise Exception('Could not parse line: %s' % line)
            summary['count'] += 1 * span
            summary['sum'] += value * span
        else:
            raise Exception('Could not parse line: %s' % line)

    return summary, mapping, write_index(summary, mapping, track)
