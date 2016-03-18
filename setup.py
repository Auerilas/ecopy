#!/usr/bin/env python

import logging
import sys
import pprint
from setuptools import setup, find_packages
from setuptools.extension import Extension

# Set up the logging environment
logging.basicConfig()
log = logging.getLogger()

# Handle the -W all flag
if 'all' in sys.warnoptions:
    log.level = logging.DEBUG

# Parse the verison from the ecopy module
with open('ecopy/__init__.py') as f:
    for line in f:
        if line.find('__version__') >= 0:
            version = line.split('=')[1].strip()
            version = version.strip('"')
            version = version.strip("'")
            continue

with open('VERSION.txt', 'w') as f:
    f.write(version)

# Use Cython if available
try:
    from Cython.Build import cythonize
except:
    log.critical(
        'Cython.Build.cythonize not found. '
        'Cython is required to build from a repo.')
    sys.exit(1)

# Use README.rst as the long description
with open('README.rst') as f:
    readme = f.read()

# Extension options
include_dirs = []
try:
    import numpy
    include_dirs.append(numpy.get_include())
except ImportError:
    log.critical('Numpy and its headers are required to run setup(). Exiting')
    sys.exit(1)

opts = dict(
    include_dirs=include_dirs,
)
log.debug('opts:\n%s', pprint.pformat(opts))

# Build extension modules 
ext_modules = cythonize([
    Extension(
        'ecopy.regression.isoFunc', ['ecopy/regression/isoFunc.pyx'], **opts),
])

# Dependencies
install_requires = [
    'numpy>=1.7',
    'scipy>=0.14',
    'matplotlib>=1.3.1',
    'pandas>=0.13',
    'patsy>=0.3.0'
]

setup_args = dict(
    name='ecopy',
    version=version,
    description='EcoPy: Ecological Data Analysis in Python',
    long_description=readme,
    url='https://github.com/Auerilas/ecopy',
    author='Nathan Lemoine',
    author_email='lemoine.nathan@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
    ],
    keywords=['ordination', 'ecology', 'multivariate data analysis'],
    ext_modules=ext_modules,
    install_requires=install_requires,
    packages=find_packages(),
)

setup(**setup_args)
