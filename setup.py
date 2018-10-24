from distutils.core import setup, Extension
from distutils.command.build import build
import os
import subprocess
import sys

class make_bwa(build):
    def run(self):
        sys.stderr.write("building BWA\n")
        os.chdir("bwa")
        subprocess.call(["make", "libbwa.a"])
        os.chdir("..")
        build.run(self)

uncalled = Extension(
    "uncalled",
     sources = ["uncalledmodule.cpp", 
                "mapper.cpp", 
                "seed_tracker.cpp", 
                "kmer_model.cpp", 
                "bwa_fmi.cpp", 
                "event_detector.cpp", 
                "range.cpp"],
     include_dirs = ["/pybind11/include", 
                     "./fast5/src", 
                     "/home-net/home-4/skovaka1@jhu.edu/anaconda3/include"],
     library_dirs = ["./bwa"],
     extra_link_args = ["-lbwa", "-lhdf5", "-lz", "-ldl"],
     extra_compile_args = ["-std=c++11"]
)

setup(name = "uncalled",
      version = "1.0",
      description = "It\"s the best around, no one\"s ever gonna get it down",
      ext_modules = [uncalled],
      cmdclass={'build': make_bwa})
