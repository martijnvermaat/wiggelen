import os
from setuptools import setup
import sys

if sys.version_info < (2, 6):
    raise Exception('Wiggelen requires Python 2.6 or higher.')

install_requires = []

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    install_requires.append('argparse')

# Python 2.6 does not include OrderedDict.
try:
    from collections import OrderedDict
except ImportError:
    install_requires.append('ordereddict')

try:
    with open('README.rst') as readme:
        long_description = readme.read()
except IOError:
    long_description = 'See https://pypi.python.org/pypi/wiggelen'

# This is quite the hack, but we don't want to import our package from here
# since that's recipe for disaster (it might have some uninstalled
# dependencies, or we might import another already installed version).
distmeta = {}
for line in open(os.path.join('wiggelen', '__init__.py')):
    try:
        field, value = (x.strip() for x in line.split('='))
    except ValueError:
        continue
    if field == '__version_info__':
        value = value.strip('[]()')
        value = '.'.join(x.strip(' \'"') for x in value.split(','))
    else:
        value = value.strip('\'"')
    distmeta[field] = value

setup(
    name='wiggelen',
    version=distmeta['__version_info__'],
    description='Working with wiggle tracks in Python',
    long_description=long_description,
    author=distmeta['__author__'],
    author_email=distmeta['__contact__'],
    url=distmeta['__homepage__'],
    license='MIT License',
    platforms=['any'],
    packages=['wiggelen'],
    install_requires=install_requires,
    entry_points = {
        'console_scripts': ['wiggelen = wiggelen.commands:main']
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
