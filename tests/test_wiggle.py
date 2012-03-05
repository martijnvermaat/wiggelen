"""
Tests for the wiggle module.
"""


import os

from nose.tools import *

import wiggelen
from wiggelen.index import INDEX_SUFFIX


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def open_(filename, mode='r'):
    """
    Open a file from the test data.
    """
    return open(os.path.join(DATA_DIR, filename), mode)


def remove_indices():
    """
    Cleanup any index files for the test data.
    """
    for file in os.listdir(DATA_DIR):
        if file.endswith(INDEX_SUFFIX):
            os.unlink(os.path.join(DATA_DIR, file))


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
