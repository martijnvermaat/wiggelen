"""
Index regions/chromosomes in wiggle tracks for random access.

Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


def index(track):
    """
    Todo: Use chaining of walker iterators to implement the use of index.
    """
    regions = {}
    while True:
        line = track.readline()
        if not line:
            break
        if 'chrom=' not in line:
            continue
        region = line.split('chrom=')[1].split()[0]
        regions[region] = track.tell() - len(line)
    return regions
