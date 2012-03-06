"""
Command line interface for working with wiggle tracks.


Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
Copyright (c) 2012 Jeroen Laros <j.f.j.laros@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""


import sys
import argparse

from .wiggle import walk, write, map_
from .index import index, write_index
from .merge import merge, mergers
from .distance import metrics, distance


def main():
    """
    Command line interface.

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

    cparser = subparsers.add_parser('scale',
        description='Scale values in a wiggle track.',
        help='scale values in a wiggle track')
    cparser.add_argument('track', metavar='TRACK',
        type=argparse.FileType('r'), help='wiggle track')
    cparser.add_argument('-f', '--factor', dest='factor', type=float,
        default=0.1, help='scaling factor to use (default: %(default)s)')

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
    dparser.add_argument('-m', dest='metric', choices=metrics, default='a',
        help='pairwise distance metric to use (default: %(default)s)')
    dparser.add_argument('tracks', metavar='TRACK', nargs='+',
        type=argparse.FileType('r'), help='wiggle track')

    args = parser.parse_args()

    if args.subcommand == 'index':
        summary, mapping, filename = index(args.track, force=True)
        if filename is None:
            parser.error('Could not write index file')

    if args.subcommand == 'sort':
        write(walk(args.track, force_index=True))

    if args.subcommand == 'scale':
        scale = lambda (r, p, v): (r, p, v * args.factor)
        write(map_(scale, walk(args.track)))

    if args.subcommand == 'merge':
        walkers = [walk(track, force_index=not args.no_indices)
                   for track in args.tracks]
        write(merge(*walkers, merger=mergers[args.merger]))

    if args.subcommand == 'distance':
        # Todo: Cleanup this code.
        distances = distance(*args.tracks, metric=metrics[args.metric])
        def name(index):
            return chr(ord('A') + index)
        try:
            sys.stdout.write(''.join('%s: %s\n' % (name(i), track.name)
                                     for i, track in enumerate(args.tracks)))
        except IOError:
            pass
        sys.stdout.write('\n   ')
        sys.stdout.write(' '.join('   %s ' % name(i)
                                  for i in range(len(args.tracks))))
        sys.stdout.write('\n')
        sys.stdout.write(name(0) + '     x\n')
        for i in range(1, len(args.tracks)):
            sys.stdout.write('%s  ' % name(i))
            for j in range(0, i):
                sys.stdout.write(' %.3f' % distances[i, j])
            sys.stdout.write('   x\n')


if __name__ == '__main__':
    main()
