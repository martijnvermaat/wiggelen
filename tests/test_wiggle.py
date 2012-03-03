"""
Tests for the wiggle module.
"""


import os

from nose.tools import *

import wiggelen


def open_(filename, mode='r'):
    dir = os.path.join(os.path.dirname(__file__), 'data')
    return open(os.path.join(dir, filename), mode)
    

class TestWiggle(object):
    """
    Tests for the wiggle module.
    """
    def test_walk(self):
        """
        Simple walk.
        """
        a = [(1, 520.0),
             (2, 536.0),
             (3, 553.0),
             (4, 568.0),
             (5, 598.0),
             (6, 616.0),
             (7, 629.0),
             (8, 649.0),
             (9, 657.0),
             (10, 676.0)]
        walker = wiggelen.walk(open_('a.mt.wig'))
        for expected, item in zip(a, walker):
            assert_equal(expected, item)
        assert_raises(StopIteration, next, walker)
