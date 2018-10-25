from distutils.core import setup, Extension
from distutils.command.build import build
import os
import subprocess
import sys

class make_bwa(build):
    def run(self):
        sys.stderr.write("building libbwa\n")
        os.chdir("bwa")
        subprocess.call(["make", "libbwa.a"])
        os.chdir("..")
        build.run(self)

uncalled = Extension(
    "uncalled",
     sources = ["src/uncalled.cpp", 
                "src/mapper.cpp", 
                "src/seed_tracker.cpp", 
                "src/kmer_model.cpp", 
                "src/bwa_fmi.cpp", 
                "src/event_detector.cpp", 
                "src/range.cpp"],
     include_dirs = ["./",
                     "./pybind11/include", 
                     "./fast5/src"],
     library_dirs = ["./bwa"],
     libraries = ["hdf5", "bwa", "z", "dl"],
     extra_compile_args = ["-std=c++11"]
)

setup(name = "uncalled",
      version = "1.0",
      description = "Rapidly maps raw nanopore signal to DNA references",
      ext_modules = [uncalled],
      cmdclass={'build': make_bwa})
