from distutils.core import setup, Extension
from distutils.command.build import build
from distutils.sysconfig import get_python_inc
import os
import subprocess
import sys

class make_bwa(build):
    def run(self):
        sys.stderr.write("building libbwa\n")
        subprocess.call(["make", 
                         "-C", "./bwa", 
                         "-f", "../src/Makefile_bwa"])
        build.run(self)

uncalled = Extension(
    "uncalled.align",
     sources = ["src/uncalled.cpp", 
                "src/mapper.cpp", 
                "src/seed_tracker.cpp", 
                "src/kmer_model.cpp", 
                "src/bwa_fmi.cpp", 
                "src/event_detector.cpp", 
                "src/range.cpp",
                "src/self_align_ref.cpp"],
     include_dirs = [get_python_inc(),
                     "./",
                     "./pybind11/include", 
                     "./fast5/include"],
     library_dirs = ["./bwa"],
     libraries = ["hdf5", "bwa", "z", "dl"],
     extra_compile_args = ["-std=c++11"]
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
      cmdclass={'build': make_bwa})
