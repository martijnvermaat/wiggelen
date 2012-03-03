import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception('Wiggelen requires Python 2.6 or higher.')

# Todo: How does this play with pip freeze requirement files?
requires = ['nose']

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

import wiggelen as distmeta

setup(
    name='Wiggelen',
    version=distmeta.__version__,
    description='Working with wiggle tracks in Python',
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license='MIT License',
    platforms=['any'],
    packages=['wiggelen'],
    requires=requires,
    entry_points = {
        'console_scripts': ['wiggelen = wiggelen.scripts:main']
        },
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        ],
    keywords='bioinformatics'
)
