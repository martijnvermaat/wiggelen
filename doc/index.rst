Wiggelen
========

Working with wiggle tracks in Python.

.. warning:: This is a work in progress, probably not yet ready for use!

.. toctree::
   :maxdepth: 2

   manual
   api


Description
-----------

Wiggelen is a Python library for working with `wiggle tracks <https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html>`_
(WIG files). It also provides a command line interface to some of its
functionality.

The main goal of Wiggelen is to provide `light-weight`_ and `unified`_ access
to wiggle tracks.

Example::

    >>> import wiggelen
    >>> for x in wiggelen.walk(open('test.wig')):
    ...     print 'chr%s:%d\t%s' % x
    chr18:34344  629.0
    chr18:34345  649.0
    chr18:34446  657.0
    chrM:308     520.0
    chrM:309     519.0


Light-weight
------------

Working with simple data should be simple. Wiggelen tries not to over-engineer
by using builtin datatypes such as tuples instead of custom objects. Sane
defaults are used throughout and things like indices are handled in the
background transparently.


Unified
-------

The central operation in Wiggelen is walking a track. Be it in fixedSteps or
variableSteps format, using any window size and step interval, walking a track
yields values one position at a time. Many operations accept walkers as input
and/or return walkers as output.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
