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
from .index import index, ordered_regions


def main():
    """
    Merge any number of wiggle tracks in various ways.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand',
        help='subcommand help')

    mparser = subparsers.add_parser('merge',
        help='merge any number of wiggle tracks in various ways')
    mparser.add_argument('-m', dest='merger', choices=mergers, default='sum',
        help='merge operation to use (default: %(default)s)')
    mparser.add_argument('tracks', metavar='TRACK', nargs='+',
        type=argparse.FileType('r'), help='wiggle track to merge')

    iparser = subparsers.add_parser('index', help='index wiggle track')
    iparser.add_argument('track', metavar='TRACK',
        type=argparse.FileType('r'), help='wiggle track to index')

    args = parser.parse_args()

    if args.subcommand == 'merge':
        #walkers = map(walk, args.tracks)
        indices = map(idx, args.tracks)
        walkers = [walk(track, regions=ordered_regions(*indices), index=idx)
                   for track, idx in zip(args.tracks, indices)]
        write(merge(*walkers, merger=mergers[args.merger]), track=sys.stdout)

    if args.subcommand == 'index':
        print index(args.track)


if __name__ == '__main__':
    main()
