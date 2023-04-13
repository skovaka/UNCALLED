from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from glob import glob
import os
import subprocess
import sys

try:
    from pybind11.setup_helpers import Pybind11Extension, ParallelCompile, naive_recompile
    ParallelCompile("NPY_NUM_BUILD_JOBS", default=4, needs_recompile=naive_recompile).install()
except:
    from setuptools import Extension as Pybind11Extension

about = {}
with open("uncalled/__about__.py") as fp:
    exec(fp.read(), about)

ROOT_DIR = os.getcwd()

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        from pybind11.setup_helpers import Pybind11Extension
        return pybind11.get_include()

uncalled = Pybind11Extension(
    "_uncalled",

    sources = glob("src/**/*.cpp", recursive=True),

    include_dirs = [
        get_pybind_include()
    ],

    libraries = ["z", "dl", "m"],

    extra_compile_args = ["-std=c++11", "-O3", "-g"],

    define_macros = [("PYBIND", None)]#, ("PYDEBUG", None)]
)

requires=[
    'pybind11>=2.6.0', 
    'pandas>=1.1.5',
    'plotly>=5.0.0',
    'dash>=2.0.0',
    'scipy>=1.5.4',
    'toml>=0.10.2',
    'ont_fast5_api',
    'pysam', 'pyslow5', 'pod5'
],

setup(
    name = about["__title__"],
    version = about["__version__"],
    description = about["__summary__"],
    author = about["__author__"],
    author_email = about["__email__"],
    url = about["__uri__"],

    classifiers=[
      "Programming Language :: Python :: 3"
    ],

    python_requires=">=3.8",

    setup_requires=['setuptools', 'pybind11>=2.6.0'],
    install_requires=requires,

    packages=find_packages(),
    include_package_data=True,
    ext_modules = [uncalled],
    #cmdclass={'build_ext': pre_build},
    entry_points = {
        'console_scripts' : ['uncalled=uncalled.__main__:main'],
    }
)
