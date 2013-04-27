Wiggelen
========

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


User documentation
------------------

New users should probably start here.

.. toctree::
   :maxdepth: 2

   install
   guide
   commands


API reference
-------------

Documentation on a specific function, class or method can be found in the API
reference.

.. toctree::
   :maxdepth: 2

   api


Additional notes
----------------

.. toctree::
   :maxdepth: 2

   development
   changelog
   copyright


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _wiggle: https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html
