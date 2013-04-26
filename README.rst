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


Documentation
-------------

The latest `documentation`_ including a user's guide and API reference is
hosted at Read The Docs.

You can also compile the documentation directly from the source code by
running ``make html`` from the ``doc/`` subdirectory. This requires `Sphinx`_
to be installed.


Source code
-----------

The Wiggelen `source code`_ is hosted on GitHub. All `unit tests`_ are run
automatically on Python 2.6, 2.7, 3.2, 3.3, and PyPy using the Travis CI
service.


Copyright
---------

Wiggelen is licensed under the MIT License, see the ``LICENSE`` file for
details. See the ``AUTHORS`` file for a list of authors.


.. _documentation: http://wiggelen.readthedocs.org/
.. _source code: https://github.com/martijnvermaat/wiggelen
.. _Sphinx: http://sphinx-doc.org/
.. _unit tests: https://travis-ci.org/martijnvermaat/wiggelen
.. _wiggle: https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html
