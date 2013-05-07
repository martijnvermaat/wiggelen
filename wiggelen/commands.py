"""
Wiggelen command line interface.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


import sys
import argparse

from .wiggle import fill, walk, write
from .index import index, write_index
from .merge import merge, mergers
from .distance import metrics, distance
from .transform import (backward_divided_difference,
                        forward_divided_difference,
                        central_divided_difference)
from . import intervals

# Python 3 compatibility.
try:
    from itertools import imap, ifilter
    map_ = imap
    filter_ = ifilter
except ImportError:
    map_ = map
    filter_ = filter

# Matplotlib only if it is installed.
try:
    from matplotlib import pyplot
except ImportError:
    pyplot = None


def log(message):
    sys.stderr.write(message + '\n')


def abort(message=None):
    if message:
        log('error: ' + message)
    sys.exit(1)


def index_track(track):
    """
    Build index for wiggle track.
    """
    # Todo: This will not rebuild the index if it already exists.
    idx, filename = index(track, force=True)
    if filename is None:
        abort('Could not write index file')


def sort_track(track):
    """
    Sort wiggle track regions alphabetically.
    """
    write(walk(track, force_index=True))


def scale_track(track, factor=0.1):
    """
    Scale values in a wiggle track.
    """
    scale = lambda (r, p, v): (r, p, v * factor)
    write(map_(scale, walk(track)))


def derivative_track(track, method='forward', step=None, auto_step=False):
    """
    Create derivative of a wiggle track.
    """
    kwargs = {'step': step}
    if method == 'central':
        derivative = central_divided_difference
    elif method == 'backward':
        derivative = backward_divided_difference
        kwargs['auto_step'] = auto_step
    else:
        derivative = forward_divided_difference
        kwargs['auto_step'] = auto_step
    write(derivative(walk(track), **kwargs))


def visualize_track(track):
    """
    Visualize a wiggle track.
    """
    # Todo: Only visualize one chromosome.
    pyplot.plot([v for _, _, v in fill(walk(track), filler=0)])
    pyplot.show()


def coverage_track(track, threshold=None):
    """
    Create coverage BED track of a wiggle track.
    """
    # Todo: Define coverage per region, like in `coverage-wiggle-to-bed` from
    #     bio-playground (https://github.com/martijnvermaat/bio-playground).
    walker = walk(track)
    if threshold is not None:
        walker = filter_(lambda (r, p, v): v >= threshold, walker)

    intervals.write(intervals.coverage(walker))


def merge_tracks(tracks, merger='sum', no_indices=False):
    """
    Merge any number of wiggle tracks in various ways.
    """
    walkers = [walk(track, force_index=not no_indices)
               for track in tracks]
    write(merge(*walkers, merger=mergers[merger]))


def distance_tracks(tracks, metric='a', threshold=None):
    """
    Calculate the distance between wiggle tracks.
    """
    # Todo: Cleanup this code.
    distances = distance(*tracks, metric=metrics[metric],
                          threshold=threshold)
    def name(index):
        return chr(ord('A') + index)
    try:
        sys.stdout.write(''.join('%s: %s\n' % (name(i), track.name)
                                 for i, track in enumerate(tracks)))
    except IOError:
        pass
    sys.stdout.write('\n   ')
    sys.stdout.write(' '.join('   %s ' % name(i)
                              for i in range(len(tracks))))
    sys.stdout.write('\n')
    sys.stdout.write(name(0) + '     x\n')
    for i in range(1, len(tracks)):
        sys.stdout.write('%s  ' % name(i))
        for j in range(0, i):
            sys.stdout.write(' %.3f' % distances[i, j])
        sys.stdout.write('   x\n')


def main():
    """
    Wiggelen command line interface.
    """
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    subparsers = parser.add_subparsers(
        title='subcommands', dest='subcommand', help='subcommand help')

    p = subparsers.add_parser(
        'index', help='build index for wiggle track',
        description=index_track.__doc__.split('\n\n')[0])
    p.set_defaults(func=index_track)
    p.add_argument(
        'track', metavar='TRACK', type=argparse.FileType('r'),
        help='wiggle track')

    p = subparsers.add_parser(
        'sort', help='sort wiggle track regions alphabetically',
        description=sort_track.__doc__.split('\n\n')[0])
    p.set_defaults(func=sort_track)
    p.add_argument(
        'track', metavar='TRACK', type=argparse.FileType('r'),
        help='wiggle track')

    p = subparsers.add_parser(
        'scale', help='scale values in a wiggle track',
        description=scale_track.__doc__.split('\n\n')[0])
    p.set_defaults(func=scale_track)
    p.add_argument(
        'track', metavar='TRACK', type=argparse.FileType('r'),
        help='wiggle track')
    p.add_argument(
        '-f', '--factor', dest='factor', type=float, default=0.1,
        help='scaling factor to use (default: %(default)s)')

    p = subparsers.add_parser(
        'derivative', help='create derivative of a wiggle track',
        description=derivative_track.__doc__.split('\n\n')[0])
    p.set_defaults(func=derivative_track)
    p.add_argument(
        'track', metavar='TRACK', type=argparse.FileType('r'),
        help='wiggle track')
    p.add_argument(
        '-m', '--method', dest='method', type=str, default='forward',
        choices=('forward', 'backward', 'central'),
        help='type of divided difference method to use (default: %(default)s)')
    p.add_argument(
        '-s', '--step', dest='step', type=int, default=None,
        help='restrict to positions that are this far apart (default: no restriction)')
    p.add_argument(
        '-a', '--auto-step', dest='auto_step', action='store_true',
        help='automatically set STEP to a value based on the first two '
        'positions in TRACK (only used if STEP is omitted, always set if '
        'METHOD is central)')

    if pyplot is not None:
        p = subparsers.add_parser(
            'visualize', help='visualize a wiggle track (requires matplotlib)',
            description=visualize_track.__doc__.split('\n\n')[0])
        p.set_defaults(func=visualize_track)
        p.add_argument(
            'track', metavar='TRACK', type=argparse.FileType('r'),
            help='wiggle track')

    p = subparsers.add_parser(
        'coverage', help='create coverage BED track of a wiggle track',
        description=coverage_track.__doc__.split('\n\n')[0])
    p.set_defaults(func=coverage_track)
    p.add_argument(
        'track', metavar='TRACK', type=argparse.FileType('r'),
        help='wiggle track')
    p.add_argument(
        '-t', '--threshold', dest='threshold', type=int, default=None,
        help='only include positions with this value or higher (default: no '
        'threshold)')

    p = subparsers.add_parser(
        'merge', help='merge any number of wiggle tracks in various ways',
        description=merge_tracks.__doc__.split('\n\n')[0])
    p.set_defaults(func=merge_tracks)
    p.add_argument(
        '-m', dest='merger', choices=mergers, default='sum',
        help='merge operation to use (default: %(default)s)')
    p.add_argument(
        '-n', '--no-indices', dest='no_indices', action='store_true',
        help='assume tracks are sorted, don\'t force building indices')
    p.add_argument(
        'tracks', metavar='TRACK', nargs='+', type=argparse.FileType('r'),
        help='wiggle track')

    # Todo: Add additional information on the metrics (using the epilog
    # argument of the subparser).
    p = subparsers.add_parser(
        'distance', help='calculate the distance between wiggle tracks',
        description=distance_tracks.__doc__.split('\n\n')[0])
    p.set_defaults(func=distance_tracks)
    p.add_argument(
        '-m', dest='metric', choices=metrics, default='a',
        help='pairwise distance metric to use (default: %(default)s)')
    p.add_argument(
        '-t', dest='threshold', type=float, default=None,
        help='threshold for noise filter (default: no noise filter)')
    p.add_argument(
        'tracks', metavar='TRACK', nargs='+', type=argparse.FileType('r'),
        help='wiggle track')

    args = parser.parse_args()

    try:
        args.func(**{k: v for k, v in vars(args).items()
                     if k not in ('func', 'subcommand')})
    except IOError as e:
        abort(str(e))


if __name__ == '__main__':
    main()
