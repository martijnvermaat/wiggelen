"""
Command line interface for working with wiggle tracks.

Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


import sys
import argparse

from . import walk, write
from .merge import merge, mergers


def main():
    """
    Merge any number of wiggle tracks in various ways.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    parser.add_argument('-m', dest='merger', choices=mergers, default='sum',
        help='merge operation to use (default: %(default)s)')
    parser.add_argument('tracks', metavar='TRACK', nargs='+',
        type=argparse.FileType('r'), help='wiggle track to merge')
    args = parser.parse_args()

    walkers = map(walk, args.tracks)
    write(merge(*walkers, merger=mergers[args.merger]), track=sys.stdout)


if __name__ == '__main__':
    main()
