"""
Tests for the distance module.
"""


import os

from nose.tools import *

from wiggelen.distance import distance
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


class TestDistance(object):
    """
    Tests for the distance module.
    """
    @classmethod
    def setup_class(cls):
        remove_indices()

    def teardown(self):
        remove_indices()

    def test_distance(self):
        """
        Simple distance test on two tracks.
        """
        track_a = open_('a.wig')
        track_b = open_('b.wig')
        assert_equal('%.3f' % distance(track_a, track_b)[1, 0],
                     '0.687')

    def test_distance_threshold(self):
        """
        Distance test with threshold on two tracks.
        """
        track_a = open_('a.wig')
        track_b = open_('b.wig')
        assert_equal('%.3f' % distance(track_a, track_b, threshold=600)[1, 0],
                     '0.925')
