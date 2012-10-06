Wiggelen
========

Working with wiggle tracks in Python.

.. warning:: This is a work in progress, probably not yet ready for use!

Wiggelen is a Python library for working with `wiggle tracks <https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html>`_
(WIG files). It also provides a command line interface to some of its
functionality.

The main goal of Wiggelen is to provide light-weight and unified access to
wiggle tracks.

Example::

    >>> import wiggelen
    >>> for x in wiggelen.walk(open('test.wig')):
    ...     print 'chr%s:%d\t%s' % x
    chr18:34344  629.0
    chr18:34345  649.0
    chr18:34446  657.0
    chrM:308     520.0
    chrM:309     519.0


Contents
--------

.. toctree::
   :maxdepth: 2

   copyright
   manual
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
