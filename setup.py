from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import os
import subprocess
import sys

__version__ = "2.2.1"

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

uncalled = Extension(
    "_uncalled",

    sources = [
       "src/pybinder.cpp",
       "src/event_profiler.cpp", 
       "src/client_sim.cpp",
       "src/fast5_reader.cpp",
       "src/mapper.cpp",
       "src/self_align_ref.cpp",
       "src/map_pool.cpp",
       "src/event_detector.cpp", 
       "src/read_buffer.cpp",
       "src/chunk.cpp",
       "src/realtime_pool.cpp",
       "src/seed_tracker.cpp", 
       "src/normalizer.cpp", 
       "src/range.cpp"
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

    define_macros = [("PYBIND", None)]
)

setup(
    name = "uncalled",
    version = __version__,
    description = "Rapidly maps raw nanopore signal to DNA references",
    author = "Sam Kovaka",
    author_email = "skovaka@gmail.com",
    url = "https://github.com/skovaka/UNCALLED",

    classifiers=[
      "Programming Language :: Python :: 3"
    ],

    python_requires=">=3.6",

    setup_requires=[
        'pybind11>=2.5.0', 
        'read-until==3.0.0'
    ],

    packages=find_packages(),
    include_package_data=True,
    scripts = ['scripts/uncalled'],
    ext_modules = [uncalled],
    cmdclass={'build_ext': pre_build},
)
