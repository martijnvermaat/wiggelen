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


Installation
------------

The Wiggelen source code is `hosted on GitHub
<https://github.com/martijnvermaat/wiggelen>`_. Supported Python versions for
running Wiggelen are 2.6, 2.7, 3.2, 3.3, and PyPy (unit tests are run
automatically on these platforms `using the Travis CI service
<https://travis-ci.org/martijnvermaat/wiggelen>`_). Wiggelen can be installed
either via the Python Package Index (PyPi) or from the source code.


Latest release via PyPi
^^^^^^^^^^^^^^^^^^^^^^^

To install the latest release via PyPi using ``pip``::

    pip install wiggelen


Development version
^^^^^^^^^^^^^^^^^^^

You can also clone and use the latest development version directly from the
GitHub repository::

    git clone https://github.com/martijnvermaat/wiggelen.git
    cd wiggelen
    python setup.py install


Documentation
-------------

The latest documentation including a user's guide and API reference is `hosted
at Read The Docs <http://wiggelen.readthedocs.org/>`_.

You can also compile the documentation directly from the source code by
running ``make html`` from the ``doc/`` subdirectory. This requires `Sphinx`_
to be installed.


Copyright
---------

Wiggelen is licensed under the MIT License, see the ``LICENSE`` file for
details. See the ``AUTHORS`` file for a list of authors.


.. _Sphinx: http://sphinx-doc.org/
.. _wiggle: https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html
