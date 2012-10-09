Wiggelen
========

Working with wiggle tracks in Python.

**Warning:** This is a work in progress, probably not yet ready for use!


Description
-----------

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


Documentation
-------------

The `latest documentation <http://wiggelen.readthedocs.org/>`_ with a manual
and API reference is hosted at Read The Docs.


Copyright
---------

Wiggelen is licensed under the MIT License, see the LICENSE file for details.
See the AUTHORS file for a list of authors.
