"""
Wiggelen, working with wiggle tracks in Python.

The `wiggle`_ (WIG) format is for storing dense, continuous genomic data such
as GC percent, probability scores, read depth, and transcriptome data.

.. _wiggle: https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html

.. todo:: Specify region(s) to use (same as filter?), possibly with a BED file.
.. todo:: Support for BigWig?

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


from .parse import ParseError
from .wiggle import ReadError, walk, zip_, fill, write


# On the event of a new release, we update the __version_info__ package
# global and set RELEASE to True.
# Before a release, a development version is denoted by a __version_info__
# ending with a 'dev' item and RELEASE is set to False.
#
# We follow a versioning scheme compatible with setuptools [1] where the
# __version_info__ variable always contains the version of the upcomming
# release (and not that of the previous release), post-fixed with a 'dev'
# item. Only in a release commit, this 'dev' item is removed (and added
# again in the next commit).
#
# [1] http://peak.telecommunity.com/DevCenter/setuptools#specifying-your-project-s-version

RELEASE = True

__version_info__ = ('0', '1', '1')
__date__ = '27 Apr 2013'


__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Martijn Vermaat, Jeroen Laros'
__contact__ = 'm.vermaat.hg@lumc.nl'
__homepage__ = 'https://github.com/martijnvermaat/wiggelen'
