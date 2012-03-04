"""
Index regions/chromosomes in wiggle tracks for random access.

Todo: Some/all of this should be private?

Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


SUFFIX = '.wix'


def write_index(index, filename):
    """
    Write the index to a file or fail silently.
    """
    try:
        with open(filename + SUFFIX, 'w') as f:
            f.write('\n'.join('%s %d' % i for i in index.items()) + '\n')
    except IOError:
        pass


def read_index(filename):
    """
    Read the index from a file or fail silently.

    Todo: Only accept if index is newer than wiggle track?
    """
    try:
        with open(filename + SUFFIX) as f:
            return dict(line.split() for line in f)
    except IOError:
        pass


def index(track):
    """
    Todo: Handle non-writable index, corrupt index, etc.
    Todo: Naming convention for index objects, other than 'index'. Perhaps use
        the same three letter name also for the index files extension? 'idx'?
    """
    filename = getattr(track, 'name', None)
    if filename is not None and filename.startswith('<'):
        filename = None

    if filename is not None:
        index = read_index(filename)
        if index is not None:
            return index

    index = {}
    while True:
        line = track.readline()
        if not line:
            break
        if 'chrom=' not in line:
            continue
        region = line.split('chrom=')[1].split()[0]
        index[region] = track.tell() - len(line)

    if filename is not None:
        write_index(index, filename)

    return index


def ordered_regions(*indices):
    """
    List of regions.
    """
    regions = set()
    for index in indices:
        regions.update(index)
    return sorted(regions)
