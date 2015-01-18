from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np

setup(
    ext_modules = cythonize("isoFunc.pyx"),
    include_dirs = [np.get_include()]
)