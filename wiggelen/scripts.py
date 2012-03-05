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
from .index import index, write_index


def main():
    """
    Merge any number of wiggle tracks in various ways.

    Todo: Organize this code based on functionality.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand',
        help='subcommand help')

    iparser = subparsers.add_parser('index',
        description='Build index for wiggle track.',
        help='build index for wiggle track')
    iparser.add_argument('track', metavar='TRACK',
        type=argparse.FileType('r'), help='wiggle track to index')

    sparser = subparsers.add_parser('sort',
        description='Sort wiggle track regions alphabetically.',
        help='sort wiggle track regions alphabetically')
    sparser.add_argument('track', metavar='TRACK',
        type=argparse.FileType('r'), help='wiggle track to sort')

    mparser = subparsers.add_parser('merge',
        description='Merge any number of wiggle tracks in various ways.',
        help='merge any number of wiggle tracks in various ways')
    mparser.add_argument('-m', dest='merger', choices=mergers, default='sum',
        help='merge operation to use (default: %(default)s)')
    mparser.add_argument('-n', '--no-indices', dest='no_indices',
        action='store_true',
        help='assume tracks are sorted, don\'t force building indices')
    mparser.add_argument('tracks', metavar='TRACK', nargs='+',
        type=argparse.FileType('r'), help='wiggle track to merge')

    args = parser.parse_args()

    if args.subcommand == 'index':
        idx = index(args.track, force=True)
        if write_index(idx, args.track) is None:
            parser.error('Could not write index file')

    if args.subcommand == 'sort':
        write(walk(args.track, force_index=True))

    if args.subcommand == 'merge':
        walkers = [walk(track, force_index=not args.no_indices)
                   for track in args.tracks]
        write(merge(*walkers, merger=mergers[args.merger]))


if __name__ == '__main__':
    main()
