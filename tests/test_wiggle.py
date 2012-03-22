"""
Tests for the wiggle module.
"""


import os
from itertools import chain

from nose.tools import *

import wiggelen
from wiggelen.index import INDEX_SUFFIX, clear_cache


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def open_(filename, mode='r'):
    """
    Open a file from the test data.
    """
    return open(os.path.join(DATA_DIR, filename), mode)


def remove_indices(keep_cache=False):
    """
    Cleanup any index files for the test data.
    """
    if not keep_cache:
        clear_cache()
    for file in os.listdir(DATA_DIR):
        if file.endswith(INDEX_SUFFIX):
            os.unlink(os.path.join(DATA_DIR, file))


def sparse(region, positions):
    """
    Create a walker for one region with given defined positions (values are
    same as position).
    """
    return ((region, p, p) for p in positions)


def filled(region, start, stop, none=[]):
    """
    Create a walker for one region with all positions between start and stop
    (inclusive). Values are same as position, except for those listed in none
    (they are None).
    """
    return ((region, p, None if p in none else p) for p in range(start, stop + 1))


class TestWiggle(object):
    """
    Tests for the wiggle module.
    """
    @classmethod
    def setup_class(cls):
        remove_indices()

    def teardown(self):
        remove_indices()

    def test_walk_single_region(self):
        """
        Walk over a track with a single region.
        """
        c = [('MT', 1, 364.0),
             ('MT', 6, 435.0),
             ('MT', 10, 485.0)]
        walker = wiggelen.walk(open_('c.wig'))
        for expected, item in zip(c, walker):
            assert_equal(expected, item)
        assert_raises(StopIteration, next, walker)

    def test_walk_multiple_regions(self):
        """
        Walk over a track with multiple regions.
        """
        values = [(2, 392.0),
                  (3, 408.0),
                  (4, 420.0),
                  (5, 452.0),
                  (7, 466.0),
                  (8, 474.0),
                  (9, 479.0)]
        b = [(r, p, v) for r in ('MT', '1', '13') for (p, v) in values]
        walker = wiggelen.walk(open_('b.wig'))
        for expected, item in zip(b, walker):
            assert_equal(expected, item)
        assert_raises(StopIteration, next, walker)

    def test_sort_multiple_regions(self):
        """
        Walk over a track with multiple regions and index.
        """
        values = [(2, 392.0),
                  (3, 408.0),
                  (4, 420.0),
                  (5, 452.0),
                  (7, 466.0),
                  (8, 474.0),
                  (9, 479.0)]
        b = [(r, p, v) for r in ('1', '13', 'MT') for (p, v) in values]
        walker = wiggelen.walk(open_('b.wig'), force_index=True)
        for expected, item in zip(b, walker):
            assert_equal(expected, item)
        assert_raises(StopIteration, next, walker)

    def test_store_index(self):
        """
        Walk over a track after the index has been made.
        """
        values = [(2, 392.0),
                  (3, 408.0),
                  (4, 420.0),
                  (5, 452.0),
                  (7, 466.0),
                  (8, 474.0),
                  (9, 479.0)]
        b = [(r, p, v) for r in ('1', '13', 'MT') for (p, v) in values]

        walker = wiggelen.walk(open_('b.wig'), force_index=True)
        for expected, item in zip(b, walker):
            assert_equal(expected, item)
        assert_raises(StopIteration, next, walker)

        walker = wiggelen.walk(open_('b.wig'))
        for expected, item in zip(b, walker):
            assert_equal(expected, item)
        assert_raises(StopIteration, next, walker)

    def test_cache_index(self):
        """
        Walk over a track after the index has been made but has been removed
        from the filesystem.
        """
        values = [(2, 392.0),
                  (3, 408.0),
                  (4, 420.0),
                  (5, 452.0),
                  (7, 466.0),
                  (8, 474.0),
                  (9, 479.0)]
        b = [(r, p, v) for r in ('1', '13', 'MT') for (p, v) in values]

        walker = wiggelen.walk(open_('b.wig'), force_index=True)
        for expected, item in zip(b, walker):
            assert_equal(expected, item)
        assert_raises(StopIteration, next, walker)

        remove_indices(keep_cache=True)

        walker = wiggelen.walk(open_('b.wig'))
        for expected, item in zip(b, walker):
            assert_equal(expected, item)
        assert_raises(StopIteration, next, walker)

    def test_walk_complex(self):
        """
        Walk over a complex track.
        """
        walker = wiggelen.walk(open_('complex.wig'))
        for _ in walker:
            pass

    def test_fill_open(self):
        """
        Test filling undefined positions.
        """
        walker = sparse('a', [3, 5, 6, 8])
        expected = list(filled('a', 3, 8, [4, 7]))
        assert_equal(list(wiggelen.fill(walker)), expected)

    def test_fill_closed(self):
        """
        Test filling undefined positions with start and stop.
        """
        walker = sparse('a', [3, 5, 6, 8])
        expected = list(filled('a', 1, 10, [1, 2, 4, 7, 9, 10]))
        assert_equal(list(wiggelen.fill(walker, regions={'a': (1, 10)})), expected)

    def test_fill_subset(self):
        """
        Test filling undefined positions on a subset.
        """
        walker = sparse('a', [1, 3, 5, 6, 8, 10])
        expected = [('a', 1, 1)] + list(filled('a', 3, 8, [4, 7])) + [('a', 10, 10)]
        assert_equal(list(wiggelen.fill(walker, regions={'a': (3, 8)})), expected)

    def test_fill_regions(self):
        """
        Test filling undefined positions over multiple regions.
        """
        a = sparse('a', [3, 5, 6, 8])
        b = sparse('b', [3, 5, 6, 8])
        c = sparse('c', [1, 3, 5, 6, 8, 10])
        walker = chain(a, b, c)
        e_a = list(sparse('a', [3, 5, 6, 8]))
        e_b = list(filled('b', 1, 10, [1, 2, 4, 7, 9, 10]))
        e_c = [('c', 1, 1)] + list(filled('c', 3, 8, [4, 7])) + [('c', 10, 10)]
        expected = list(chain(e_a, e_b, e_c))
        assert_equal(list(wiggelen.fill(walker, regions={'b': (1, 10), 'c': (3, 8)})), expected)
