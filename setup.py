from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext
from distutils.sysconfig import get_python_inc
import os
import subprocess
import sys

class make_bwa(build_ext):
    def run(self):
        sys.stderr.write("building libbwa\n")
        subprocess.call(["make", 
                         "-C", "./submods/bwa", 
                         "-f", "../../src/Makefile_bwa"])
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
                     "./submods/pybind11/include", 
                     "./submods/fast5/include",
                     "./submods/pdqsort",
                     "./submods/toml11"],
     #library_dirs = ["./submods/bwa ./submods/bwa/libbwa.a", "./submods/hdf5/lib ./submods/hdf5/lib/libhdf5.a"],
     library_dirs = ["./submods/bwa", "./submods/hdf5/lib"],
     libraries = ["bwa", "hdf5", "z", "dl", "m"],
     extra_compile_args = ["-std=c++11", "-O3"],
     define_macros = [("PYBIND", None)]
)

setup(name = "uncalled",
      version = "1.3",
      description = "Rapidly maps raw nanopore signal to DNA references",
      author = "Sam Kovaka",
      author_email = "skovaka@gmail.com",
      url = "https://github.com/skovaka/UNCALLED",
      packages=['uncalled'],
      package_dir = {'': 'src', 'uncalled': 'src/uncalled'},
      py_modules = ['uncalled.params', 'uncalled.minknow_client', 'uncalled.index', 'uncalled.pafstats'],
      ext_modules = [uncalled],
      package_data = {'uncalled': ['models/*', 'conf/*']},
      scripts = ['scripts/uncalled'],
      cmdclass={'build_ext': make_bwa})
