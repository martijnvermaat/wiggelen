"""
Index regions/chromosomes in wiggle tracks for random access.

Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


import sys


INDEX_SUFFIX = '.idx'


def index_filename(track=sys.stdin):
    filename = getattr(track, 'name', None)
    if filename is not None and not filename.startswith('<'):
        return filename + INDEX_SUFFIX


def write_index(idx, track=sys.stdout):
    """
    Write the index to a file or fail silently.
    """
    filename = index_filename(track)

    if filename is not None:
        try:
            with open(filename, 'w') as f:
                f.write('\n'.join('%s %d' % i for i in idx.items()) + '\n')
                return filename
        except IOError:
            pass


def read_index(track=sys.stdin):
    """
    Read the index from a file or fail silently and return None.

    Todo: Only accept if index is newer than wiggle track?
    """
    filename = index_filename(track)

    if filename is not None:
        try:
            with open(filename) as f:
                return dict((r, int(p)) for r, p in (l.split() for l in f))
        except IOError:
            pass


def index(track=sys.stdin, force=False):
    """
    Return index of region positions in track.

    @arg track: Wiggle track (default: standard input).
    @type track: file
    @kwarg force: Force creating an index if it does not yet exist
        (default: False).
    @type force: bool

    Todo: Handle non-writable index, corrupt index, etc.
    Todo: Also including the end positions would make it possible to do random
        jumps inside a region with some educated guessing. Perfect hits would
        not be possible, since the length of the lines are variable.
    """
    idx = read_index(track)

    if idx is not None or not force:
        return idx

    idx = {}
    while True:
        line = track.readline()
        if not line:
            break
        if 'chrom=' not in line:
            continue
        region = line.split('chrom=')[1].split()[0]
        idx[region] = track.tell() - len(line)

    write_index(idx, track)
    return idx
