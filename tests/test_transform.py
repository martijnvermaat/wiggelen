"""
Tests for the transform module.
"""


import os
from itertools import chain

from nose.tools import *

from wiggelen.transform import (forward_divided_difference,
                                backward_divided_difference,
                                central_divided_difference)


class TestTransform(object):
    """
    Tests for the transform module.

    Todo: More tests for divided difference functions.
    """
    def test_forward_divided_difference(self):
        """
        Simple forward divided difference.
        """
        orig = [(1, 5), (2, 4), (3, 4), (4, 4), (5, 5), (6, 4), (7, 3), (8, 1), (9, 5), (10, 6)]
        expected = [(1, -1), (2, 0), (3, 0), (4, 1), (5, -1), (6, -1), (7, -2), (8, 4), (9, 1)]
        walker = (('a', p, v) for p, v in orig)
        assert_equal([(p, v) for _, p, v in forward_divided_difference(walker)], expected)

    def test_backward_divided_difference(self):
        """
        Simple backward divided difference.
        """
        orig = [(1, 5), (2, 4), (3, 4), (4, 4), (5, 5), (6, 4), (7, 3), (8, 1), (9, 5), (10, 6)]
        expected = [(2, -1), (3, 0), (4, 0), (5, 1), (6, -1), (7, -1), (8, -2), (9, 4), (10, 1)]
        walker = (('a', p, v) for p, v in orig)
        assert_equal([(p, v) for _, p, v in backward_divided_difference(walker)], expected)

    def test_central_divided_difference(self):
        """
        Simple central divided difference.
        """
        orig = [(1, 5), (2, 4), (3, 4), (4, 4), (5, 5), (6, 4), (7, 3), (8, 1), (9, 5), (10, 6)]
        expected = [(2, -0.5), (3, 0), (4, 0.5), (5, 0), (6, -1), (7, -1.5), (8, 1), (9, 2.5)]
        walker = (('a', p, v) for p, v in orig)
        assert_equal([(p, v) for _, p, v in central_divided_difference(walker)], expected)
