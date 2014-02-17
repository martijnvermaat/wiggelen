"""
Wiggelen, working with wiggle tracks in Python.

The `wiggle`_ (WIG) format is for storing dense, continuous genomic data such
as GC percent, probability scores, read depth, and transcriptome data.

.. _wiggle: https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE file.
"""


from .parse import ParseError
from .wiggle import ReadError, walk, zip_, fill, write


# We follow a versioning scheme compatible with setuptools [1] where the
# package version is always that of the upcoming release (and not that of the
# previous release), post-fixed with ``.dev``. Only in a release commit, the
# ``.dev`` is removed (and added again in the next commit).
#
# Note that this scheme is not 100% compatible with SemVer [2] which would
# require ``-dev`` instead of ``.dev``.
#
# [1] http://peak.telecommunity.com/DevCenter/setuptools#specifying-your-project-s-version
# [2] http://semver.org/


__version_info__ = ('0', '4', '0')
__date__ = '17 Feb 2014'


__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Martijn Vermaat, Jeroen Laros'
__contact__ = 'm.vermaat.hg@lumc.nl'
__homepage__ = 'https://github.com/martijnvermaat/wiggelen'
