"""
Wiggelen command line interface.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


from __future__ import division

import argparse
import sys

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

# Plotting only if Matplotlib is installed.
try:
    from matplotlib import pyplot
    from .plot import plot
except ImportError:
    matplotlib = plot = None


def log(message):
    sys.stderr.write(message + '\n')


def abort(message=None):
    if message:
        log('error: ' + message)
    sys.exit(1)


def read_regions(regions):
    """
    Ad-hoc parsing of regions from a BED track.
    """
    result = {}
    for line in regions:
        if line.startswith('track'):
            continue
        region, start, stop = line.strip().split('\t')[:3]
        result[region] = int(start) + 1, int(stop)
    return result


def index_track(track):
    """
    Build index for wiggle track.
    """
    # Todo: This will not rebuild the index if it already exists.
    idx, filename = index(track, force=True)
    if filename is None:
        abort('Could not write index file')


def sort_track(track, name=None, description=None):
    """
    Sort wiggle track regions alphabetically.
    """
    if name is None and hasattr(track, 'name'):
        name = 'Sorted %s' % track.name

    write(walk(track, force_index=True), name=name, description=description)


def scale_track(track, factor=0.1, name=None, description=None):
    """
    Scale values in a wiggle track.
    """
    if name is None and hasattr(track, 'name'):
        name = 'Scaled %s' % track.name

    scale = lambda (r, p, v): (r, p, v * factor)
    write(map_(scale, walk(track)), name=name, description=description)


def fill_track(track, genome=None, filler='0', only_edges=False, name=None,
               description=None):
    """
    Fill in undefined positions in a wiggle track.
    """
    if name is None and hasattr(track, 'name'):
        name = 'Filled %s' % track.name

    if genome is not None:
        genome = read_regions(genome)

    try:
        filler = float(filler) if '.' in filler else int(filler)
    except ValueError:
        abort('Could not parse filler value: %s' % filler)

    write(fill(walk(track),
               regions=genome, filler=filler, only_edges=only_edges),
          name=name, description=description)


def derivative_track(track, method='forward', step=None, auto_step=False,
                     name=None, description=None):
    """
    Create derivative of a wiggle track.
    """
    if name is None and hasattr(track, 'name'):
        name = 'Derivative of %s' % track.name

    kwargs = {'step': step}
    if method == 'central':
        derivative = central_divided_difference
    elif method == 'backward':
        derivative = backward_divided_difference
        kwargs['auto_step'] = auto_step
    else:
        derivative = forward_divided_difference
        kwargs['auto_step'] = auto_step
    write(derivative(walk(track), **kwargs), name=name,
          description=description)


def plot_tracks(tracks, regions=None, genome=None, order_by='region',
                average_threshold=None, columns=None, pdf=None):
    """
    Visualize wiggle tracks in a plot.
    """
    def track_name(track, i):
        return track.name if hasattr(track, 'name') else 'track %i' % i

    # Prepend track name if we have more than one track.
    def rename(region, track, i):
        if len(tracks) > 1:
            return track_name(track, i) + ', ' + region
        return region

    if genome is not None:
        # Read genome from BED track and filter by specified regions. If we
        # have more than one track, have a region copy per track.
        genome = {rename(r, track, i): v
                  for r, v in read_regions(genome).items()
                  for i, track in enumerate(tracks)
                  if regions is None or r in regions}

    # Filter by specified regions.
    def filtered(walker):
        if regions is None:
            return walker
        return filter_(lambda (r, p, v): r in regions, walker)

    # Have all tracks concatenated in one walker.
    walker = ((rename(r, track, i), p, v)
              for i, track in enumerate(tracks)
              for r, p, v in filtered(walk(track)))

    fig, axes, rows, columns = plot(
        walker, regions=genome, order_by=order_by,
        average_threshold=average_threshold, columns=columns)

    if pdf:
        fig.set_size_inches(6 * columns, 3 * rows)
        fig.savefig(pdf, format='pdf')
    else:
        pyplot.show()


def coverage_track(track, threshold=None, name=None, description=None):
    """
    Create coverage BED track of a wiggle track.
    """
    if name is None and hasattr(track, 'name'):
        name = 'Coverage of %s' % track.name

    # Todo: Define coverage per region, like in `coverage-wiggle-to-bed` from
    #     bio-playground (https://github.com/martijnvermaat/bio-playground).
    walker = walk(track)
    if threshold is not None:
        walker = filter_(lambda (r, p, v): v >= threshold, walker)

    intervals.write(intervals.coverage(walker), name=name,
                    description=description)


def merge_tracks(tracks, merger='sum', no_indices=False, name=None,
                 description=None):
    """
    Merge any number of wiggle tracks in various ways.
    """
    if name is None and all(hasattr(track, 'name') for track in tracks):
        name = 'Merge of %s' % ', '.join(track.name for track in tracks)

    walkers = [walk(track, force_index=not no_indices)
               for track in tracks]
    write(merge(*walkers, merger=mergers[merger]), name=name,
          description=description)


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
    p.add_argument(
        '-n', '--name', dest='name', type=str,
        help='name to use for result track, displayed to the left of the '
        'track in the UCSC Genome Browser (default: Sorted TRACK)')
    p.add_argument(
        '-d', '--description', dest='description', type=str,
        help='description to use for result track, displayed as center label '
        'in the UCSC Genome Browser (default: no description)')

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
    p.add_argument(
        '-n', '--name', dest='name', type=str,
        help='name to use for result track, displayed to the left of the '
        'track in the UCSC Genome Browser (default: Sorted TRACK)')
    p.add_argument(
        '-d', '--description', dest='description', type=str,
        help='description to use for result track, displayed as center label '
        'in the UCSC Genome Browser (default: no description)')

    p = subparsers.add_parser(
        'fill', help='fill undefined positions in a wiggle track',
        description=fill_track.__doc__.split('\n\n')[0],
        epilog='Note that the resulting track may be very large if '
        '--only-edges is not specified.')
    p.set_defaults(func=fill_track)
    p.add_argument(
        'track', metavar='TRACK', type=argparse.FileType('r'),
        help='wiggle track')
    p.add_argument(
        '-g', '--genome', dest='genome', type=argparse.FileType('r'),
        help='regions in BED format to fill (if not specified, first and '
        'last defined positions in the track define the regions)')
    p.add_argument(
        '-f', '--filler', dest='filler', default='0',
        help='value to use for undefined positions (default: %(default)s)')
    p.add_argument(
        '-e', '--only-edges', dest='only_edges', action='store_true',
        help='only fill the first and last of continuously undefined '
        'positions')
    p.add_argument(
        '-n', '--name', dest='name', type=str,
        help='name to use for result track, displayed to the left of the '
        'track in the UCSC Genome Browser (default: Sorted TRACK)')
    p.add_argument(
        '-d', '--description', dest='description', type=str,
        help='description to use for result track, displayed as center label '
        'in the UCSC Genome Browser (default: no description)')

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
    p.add_argument(
        '-n', '--name', dest='name', type=str,
        help='name to use for result track, displayed to the left of the '
        'track in the UCSC Genome Browser (default: Sorted TRACK)')
    p.add_argument(
        '-d', '--description', dest='description', type=str,
        help='description to use for result track, displayed as center label '
        'in the UCSC Genome Browser (default: no description)')

    if plot is not None:
        p = subparsers.add_parser(
            'plot', description=plot_tracks.__doc__.split('\n\n')[0],
            help='visualize wiggle tracks in a plot (requires matplotlib)')
        p.set_defaults(func=plot_tracks)
        p.add_argument(
            'tracks', metavar='TRACK', type=argparse.FileType('r'), nargs='+',
            help='wiggle track')
        p.add_argument(
            '-r', '--regions', dest='regions', type=str, default=None, nargs='+',
            help='plot only these regions (default: plot all regions with data)')
        p.add_argument(
            '-g', '--genome', dest='genome', type=argparse.FileType('r'),
            help='regions with start and end positions in BED format')
        p.add_argument(
            '--order-by', dest='order_by', type=str, default='region',
            choices=['region', 'average', 'original'],
            help='when plotting multiple tracks and/or regions, order them '
            'by this key (default: region)')
        p.add_argument(
            '-t', '--avg-threshold', dest='average_threshold', type=float,
            default=None, help='Plot regions with average value above or '
            'equal to this threshold (default: plot regions that have data)')
        p.add_argument(
            '-c', '--columns', dest='columns', type=int, default=None,
            help='when plotting multiple tracks and/or regions, use this many columns (default: automatically chosen)')
        p.add_argument(
            '-o', '--output', dest='pdf', type=argparse.FileType('wb'),
            default=None, help='output PDF file')

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
    p.add_argument(
        '-n', '--name', dest='name', type=str,
        help='name to use for result track, displayed to the left of the '
        'track in the UCSC Genome Browser (default: Sorted TRACK)')
    p.add_argument(
        '-d', '--description', dest='description', type=str,
        help='description to use for result track, displayed as center label '
        'in the UCSC Genome Browser (default: no description)')

    p = subparsers.add_parser(
        'merge', help='merge any number of wiggle tracks in various ways',
        description=merge_tracks.__doc__.split('\n\n')[0])
    p.set_defaults(func=merge_tracks)
    p.add_argument(
        '-m', dest='merger', choices=mergers, default='sum',
        help='merge operation to use (default: %(default)s)')
    p.add_argument(
        '-x', '--no-indices', dest='no_indices', action='store_true',
        help='assume tracks are sorted, don\'t force building indices')
    p.add_argument(
        'tracks', metavar='TRACK', nargs='+', type=argparse.FileType('r'),
        help='wiggle track')
    p.add_argument(
        '-n', '--name', dest='name', type=str,
        help='name to use for result track, displayed to the left of the '
        'track in the UCSC Genome Browser (default: Sorted TRACK)')
    p.add_argument(
        '-d', '--description', dest='description', type=str,
        help='description to use for result track, displayed as center label '
        'in the UCSC Genome Browser (default: no description)')

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
