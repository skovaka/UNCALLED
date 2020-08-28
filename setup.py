from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import subprocess
import sys

__version__ = "2.1"

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        return pybind11.get_include()

class make_libs(build_ext):
    def run(self):
        sys.stderr.write("building libbwa\n")

        subprocess.call([
            "make", 
             "-C", "./submods/bwa", 
             "-f", "../../src/Makefile_bwa"
        ])

        wd = os.getcwd()

        os.chdir("submods/hdf5")

        subprocess.call([
            "./configure", 
                "--enable-threadsafe", 
                "--disable-hl",
                "--prefix=`pwd`" ,
                "--enable-shared=no",
                "--with-pic=yes"
        ])

        subprocess.call(["make"])
        subprocess.call(["make", "install"])

        os.chdir(wd)

        build_ext.run(self)

uncalled = Extension(
    "uncalled.mapping",
     sources = [
                "src/uncalled.cpp",
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
                "src/range.cpp"],
     include_dirs = ["./submods", #TODO: consistent incl paths
                     "./submods/hdf5/include", 
                     "./submods/fast5/include",
                     "./submods/pdqsort",
                     "./submods/toml11",
                     get_pybind_include()],
     #library_dirs = ["./submods/bwa ./submods/bwa/libbwa.a", "./submods/hdf5/lib ./submods/hdf5/lib/libhdf5.a"],
     library_dirs = ["./submods/bwa", "./submods/hdf5/lib"],
     libraries = ["bwa", "hdf5", "z", "dl", "m"],
     extra_compile_args = ["-std=c++11", "-O3"],
     define_macros = [("PYBIND", None)]
)

setup(name = "uncalled",
      version = __version__,
      description = "Rapidly maps raw nanopore signal to DNA references",
      author = "Sam Kovaka",
      author_email = "skovaka@gmail.com",
      url = "https://github.com/skovaka/UNCALLED",
      install_requires=['pybind11>=2.5.0', 'numpy>=1.17.4'],
      packages=['uncalled'],
      package_dir = {'': 'src', 'uncalled': 'src/uncalled'},
      py_modules = ['uncalled.params', 'uncalled.minknow_client', 'uncalled.index', 'uncalled.pafstats'],
      ext_modules = [uncalled],
      package_data = {'uncalled': ['models/*', 'conf/*']},
      scripts = ['scripts/uncalled'],
      classifiers=[
          "Programming Language :: Python :: 3"
      ],
      cmdclass={'build_ext': make_libs},
      python_requires=">=3.6")
