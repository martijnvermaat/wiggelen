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
from .index import index, write_index
from .merge import merge, mergers
from .distance import distance


def main():
    """
    Merge any number of wiggle tracks in various ways.

    .. todo:: Organize this code based on functionality.
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
        type=argparse.FileType('r'), help='wiggle track')

    sparser = subparsers.add_parser('sort',
        description='Sort wiggle track regions alphabetically.',
        help='sort wiggle track regions alphabetically')
    sparser.add_argument('track', metavar='TRACK',
        type=argparse.FileType('r'), help='wiggle track')

    mparser = subparsers.add_parser('merge',
        description='Merge any number of wiggle tracks in various ways.',
        help='merge any number of wiggle tracks in various ways')
    mparser.add_argument('-m', dest='merger', choices=mergers, default='sum',
        help='merge operation to use (default: %(default)s)')
    mparser.add_argument('-n', '--no-indices', dest='no_indices',
        action='store_true',
        help='assume tracks are sorted, don\'t force building indices')
    mparser.add_argument('tracks', metavar='TRACK', nargs='+',
        type=argparse.FileType('r'), help='wiggle track')

    dparser = subparsers.add_parser('distance',
        description='Calculate the distances between wiggle tracks.',
        help='calculate the distances between wiggle tracks')
    dparser.add_argument('-n', '--no-indices', dest='no_indices',
        action='store_true',
        help='assume tracks are sorted, don\'t force building indices')
    dparser.add_argument('tracks', metavar='TRACK', nargs='+',
        type=argparse.FileType('r'), help='wiggle track')

    args = parser.parse_args()

    if args.subcommand == 'index':
        summary, mapping = index(args.track, force=True)
        # Note that this writes the index twice if it does not yet exist.
        if write_index(summary, mapping, args.track) is None:
            parser.error('Could not write index file')

    if args.subcommand == 'sort':
        write(walk(args.track, force_index=True))

    if args.subcommand == 'merge':
        walkers = [walk(track, force_index=not args.no_indices)
                   for track in args.tracks]
        write(merge(*walkers, merger=mergers[args.merger]))

    if args.subcommand == 'distance':
        distances = distance(*args.tracks)
        names = 'ABCDEFGH'
        sys.stdout.write('   ' + ' '.join('  %s ' % n for n in names[:len(args.tracks)]) + '\n')
        sys.stdout.write(names[0] + '    x\n')
        for i in range(1, len(args.tracks)):
            sys.stdout.write('%s  ' % names[i])
            for j in range(0, i):
                sys.stdout.write(' %.2f' % distances[i, j])
            sys.stdout.write('  x\n')


if __name__ == '__main__':
    main()
