Wiggelen
========

Working with wiggle tracks in Python.

Wiggelen is a Python library for working with `wiggle`_ tracks (WIG files). It
also provides a command line interface to some of its functionality.

The main goal of Wiggelen is to provide light-weight and unified access to
wiggle tracks.

::

    >>> import wiggelen
    >>> for x in wiggelen.walk(open('test.wig')):
    ...     print 'chr%s:%d\t%s' % x
    ...
    chr18:34344  629.0
    chr18:34345  649.0
    chr18:34446  657.0
    chrM:308     520.0
    chrM:309     519.0


Contents
--------

.. toctree::
   :maxdepth: 2

   intro
   guide
   commands
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _wiggle: https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html
