"""
Tests for the intervals module.
"""


from nose.tools import *

from wiggelen.intervals import coverage


class TestIntervals(object):
    """
    Tests for the intervals module.
    """
    def test_coverage(self):
        """
        Simple interval coverage on one region.
        """
        orig = [(1, 5), (2, 4), (3, 4), (4, 4), (5, 5), (7, 3), (8, 1), (9, 5), (10, 6)]
        expected = [(1, 5), (7, 10)]
        walker = (('a', p, v) for p, v in orig)
        assert_equal([(b, e) for _, b, e in coverage(walker)], expected)

    def test_coverage_regions(self):
        """
        Interval coverage on two regions.
        """
        orig = [('a', 1, 5), ('a', 2, 4), ('a', 4, 4), ('a', 5, 5),
                ('b', 6, 5), ('b', 7, 4), ('b', 8, 4), ('b', 12, 4)]
        expected = [('a', 1, 2), ('a', 4, 5), ('b', 6, 8), ('b', 12, 12)]
        assert_equal(list(coverage(orig)), expected)

    def test_coverage_empty(self):
        """
        Interval coverage on empty walker.
        """
        orig = []
        expected = []
        assert_equal(list(coverage(orig)), expected)
