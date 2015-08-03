from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

DESCRIPTION = 'EcoPy: Ecological Data Analysis in Python'
LONG_DESCRIPTION = """\
EcoPy is a module for multivariate data analysis in Python. It is built on numpy, pandas, and in some instances, scipy.

Some features of EcoPy are:
- Principle components analysis
- Correspondance analysis
- Redundancy analysis
- ANOSIM/SIMPER
- Numerous other techniques
"""

def check_dependencies():
    install_requires = []

    try:
        import numpy
    except ImportError:
        install_requires.append('numpy')
    try:
        import scipy
    except ImportError:
        install_requires.append('scipy')
    try:
        import matplotlib
    except ImportError:
        install_requires.append('matplotlib')
    try:
        import pandas
    except ImportError:
        install_requires.append('pandas')

    return install_requires

install_requires = check_dependencies()

setup(
    name='ecopy',
    version= '0.0.7a1',
    description = DESCRIPTION,
    long_description = LONG_DESCRIPTION,
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
    packages=['ecopy', 'ecopy.base_funcs', 'ecopy.diversity', 'ecopy.matrix_comp', 'ecopy.ordination', 'ecopy.regression'],
    install_requires=install_requires,
)