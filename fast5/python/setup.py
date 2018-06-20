#!/usr/bin/env python

#
# Part of: https://github.com/mateidavid/fast5
#
# (c) 2017: Matei David, Ontario Institute for Cancer Research
# MIT License
#

import os
import sys

from setuptools import setup, Extension

use_cython = True #os.environ.get('USE_CYTHON', '') != ''

# check HDF5 include and lib dirs
hdf5_dir = os.environ.get('HDF5_DIR', '/usr')
hdf5_include_dir = os.environ.get('HDF5_INCLUDE_DIR', os.path.join(hdf5_dir, 'include'))
hdf5_lib_dir = os.environ.get('HDF5_LIB_DIR', os.path.join(hdf5_dir, 'lib'))
hdf5_lib = os.environ.get('HDF5_LIB', 'hdf5')
if not os.path.isfile(os.path.join(hdf5_include_dir, 'H5pubconf.h')):
    sys.exit(hdf5_include_dir + ': could not find HDF5 header files; use HDF5_DIR or HDF5_INCLUDE_DIR')
if (not os.path.isfile(os.path.join(hdf5_lib_dir, 'lib' + hdf5_lib + '.so'))
    and not os.path.isfile(os.path.join(hdf5_lib_dir, 'lib' + hdf5_lib + '.a'))):
    sys.exit(hdf5_lib_dir + ': could not find HDF5 library file; use HDF5_DIR or HDF5_LIB_DIR/HDF5_LIB')

fast5_dir = os.environ.get('FAST5_DIR', '..')
fast5_src_dir = os.path.join(fast5_dir, 'src')
fast5_version = open(os.path.join(fast5_dir, 'VERSION')).readline().strip()

extra_compile_args = [
    '-std=c++11',
    '-Wall', '-Wextra', '-Wpedantic',
]
# don't indiscriminately add /usr/include to work around bug:
# https://lists.fedoraproject.org/archives/list/devel@lists.fedoraproject.org/thread/Q5SWCUUMWQ4EMS7CU2CBOZHV3WZYOOTT/
for d in [hdf5_include_dir]:
    if d != '/usr/include':
        extra_compile_args += ['-isystem', d]
#extra_compile_args += ['-O0', '-g3', '-ggdb', '-fno-eliminate-unused-debug-types', '-v']

extra_link_args = []
#extra_link_args += ['-v']

#if sys.platform == 'darwin':
#    extra_compile_args.append('-mmacosx-version-min=10.7')

extensions = [
    Extension(
        'fast5',
        language='c++',
        sources=['fast5/fast5.' + ['cpp', 'pyx'][use_cython]],
        include_dirs=[fast5_src_dir],
        library_dirs=[hdf5_lib_dir],
        runtime_library_dirs=[hdf5_lib_dir],
        libraries=[hdf5_lib],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    ),
]

if use_cython:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name='fast5',
    description='Fast5 file interface.',
    version=fast5_version,
    author='Matei David, Ontario Institute for Cancer Research',
    author_email='matei.david at oicr.on.ca',
    license='MIT',
    url='https://github.com/mateidavid/fast5',
    ext_modules=extensions,
    scripts=[
        os.path.join('bin', 'f5ls'),
        os.path.join('bin', 'f5pack'),
    ],
)
