from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from glob import glob
import os
import subprocess
import sys

try:
    from pybind11.setup_helpers import Pybind11Extension, ParallelCompile, naive_recompile
    ParallelCompile("NPY_NUM_BUILD_JOBS", needs_recompile=naive_recompile).install()
except:
    from setuptools import Extension as Pybind11Extension

about = {}
with open("uncalled/__about__.py") as fp:
    exec(fp.read(), about)

ROOT_DIR = os.getcwd()

SUBMOD_DIR = os.path.join(ROOT_DIR, "submods")
SUBMODS = [
    "bwa", 
    "fast5", 
    "hdf5", 
    "pdqsort", 
    "toml11"
]

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        from pybind11.setup_helpers import Pybind11Extension
        return pybind11.get_include()

class pre_build(build_ext):
    def run(self):

        submods_loaded = True
        for submod in SUBMODS:
            if not os.path.exists(os.path.join(SUBMOD_DIR, submod)):
                submods_loaded = False
                break

        if not submods_loaded:
            sys.stderr.write("Downloading submodules\n")
            subprocess.check_call([
                "git", "submodule", "update", "--init"
            ])
        else:
            sys.stderr.write("All submodules present\n")

        if os.path.exists("./submods/bwa/libbwa.a"):
            sys.stderr.write("Found libbwa.a\n")
        else:
            sys.stderr.write("building libbwa\n")

            subprocess.check_call([
                "make", 
                 "-C", "./submods/bwa", 
                 "-f", "../../src/Makefile_bwa"
            ])

        if os.path.exists("./submods/hdf5/lib/libhdf5.a"):
            sys.stderr.write("Found libhdf5.a\n")
        else:

            hdf5_dir = os.path.join(ROOT_DIR, "submods/hdf5")

            os.chdir(hdf5_dir)

            subprocess.check_call([
                "./configure", 
                    "--enable-threadsafe", 
                    "--disable-hl",
                    "--prefix", hdf5_dir,
                    "--enable-shared=no",
                    "--with-pic=yes"
            ])

            subprocess.check_call(["make"])
            subprocess.check_call(["make", "install"])

            os.chdir(ROOT_DIR)

        build_ext.run(self)

uncalled = Pybind11Extension(
    "_uncalled",

    sources = [ #glob("src/**/*.cpp", recursive=True),
       "src/pore_model.cpp",
       "src/dtw.cpp",
       "src/pybinder.cpp",
       "src/config.cpp",
       "src/event_profiler.cpp", 
       "src/simulator.cpp",
       "src/fast5_reader.cpp",
       "src/mapper.cpp",
       "src/self_align_ref.cpp",
       "src/map_pool.cpp",
       "src/map_pool_ord.cpp",
       "src/event_detector.cpp", 
       "src/read_buffer.cpp",
       "src/realtime_pool.cpp",
       "src/seed_tracker.cpp", 
       "src/normalizer.cpp", 
       "src/paf.cpp", 
       "src/range.cpp",
       "src/dataframe.cpp"
    ],

    include_dirs = [
        "./submods",
        "./submods/hdf5/include", 
        "./submods/fast5/include",
        "./submods/pdqsort",
        "./submods/toml11",
        get_pybind_include()
    ],

    library_dirs = [
        "./submods/bwa"
    ],
    
    extra_objects = [                 
        "./submods/hdf5/lib/libhdf5.a"
    ],                                

    libraries = ["bwa", "z", "dl", "m"],

    extra_compile_args = ["-std=c++11", "-O3"],

    define_macros = [("PYBIND", None)]#, ("PYDEBUG", None)]
)

requires=[
    'pybind11>=2.6.0', 
    #'read-until==3.0.0',
    'pandas>=1.1.5',
    'plotly>=5.0.0',
    'dash>=2.0.0',
    'scipy>=1.5.4',
    'toml>=0.10.2',
    'ont_fast5_api',
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

    setup_requires=['pybind11>=2.6.0'],
    install_requires=requires,

    packages=find_packages(),
    include_package_data=True,
    ext_modules = [uncalled],
    cmdclass={'build_ext': pre_build},
    entry_points = {
        'console_scripts' : ['uncalled=uncalled.__main__:main'],
    }
)
