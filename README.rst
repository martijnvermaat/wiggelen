Wiggelen
========

Wiggelen is a Python library for working with `wiggle`_ tracks (WIG files). It
also provides a command line interface to some of its functionality.

The main goal of Wiggelen is to provide light-weight and unified access to
wiggle tracks.

.. code-block:: pycon

    >>> import wiggelen
    >>> for x in wiggelen.walk(open('test.wig')):
    ...     print 'chr%s:%d\t%s' % x
    ...
    chr18:34344  629.0
    chr18:34345  649.0
    chr18:34446  657.0
    chrM:308     520.0
    chrM:309     519.0

To install the latest release via PyPI using pip::

    pip install wiggelen


Documentation
-------------

The `latest documentation <http://wiggelen.readthedocs.org/>`_ with user guide
and API reference is hosted at Read The Docs.

You can also compile the documentation directly from the source code by
running ``make html`` from the ``doc/`` subdirectory. This requires `Sphinx`_
to be installed.


Copyright
---------

Wiggelen is licensed under the MIT License, see the LICENSE file for
details. See the AUTHORS file for a list of authors.


.. _Sphinx: http://sphinx-doc.org/
.. _wiggle: https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html
