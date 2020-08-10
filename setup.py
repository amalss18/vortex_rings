from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Velocity field',
    ext_modules=cythonize("fasten.pyx"),
    zip_safe=False,
)