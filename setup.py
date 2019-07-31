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
                         "-C", "./bwa", 
                         "-f", "../src/Makefile_bwa"])
        build_ext.run(self)

uncalled = Extension(
    "uncalled.mapping",
     sources = ["src/mapper.cpp",
                "src/uncalled.cpp",
                "src/read_buffer.cpp",
                "src/fast5_pool.cpp",
                "src/params.cpp",
                "src/simulator.cpp",
                "src/bwa_fmi.cpp", 
                "src/chunk.cpp",
                "src/chunk_pool.cpp",
                "src/seed_tracker.cpp", 
                "src/normalizer.cpp", 
                "src/kmer_model.cpp", 
                "src/event_detector.cpp", 
                "src/range.cpp",
                "src/self_align_ref.cpp"],
     include_dirs = ["./",
                     "./pybind11/include", 
                     "./fast5/include",
                     "./pdqsort"],
     library_dirs = ["./bwa"],
     libraries = ["hdf5", "bwa", "z", "dl"],
     extra_compile_args = ["-std=c++11", "-O3"]
)

setup(name = "uncalled",
      version = "0.1",
      description = "Rapidly maps raw nanopore signal to DNA references",
      author = "Sam Kovaka",
      author_email = "skovaka@gmail.com",
      url = "https://github.com/skovaka/UNCALLED",
      packages=['uncalled'],
      package_dir = {'': 'src', 'uncalled': 'src/uncalled'},
      py_modules = ['uncalled.index'],
      ext_modules = [uncalled],
      package_data = {'uncalled': ['models/*']},
      scripts = ['scripts/uncalled'],
      cmdclass={'build_ext': make_bwa})
