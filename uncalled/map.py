#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import division
import sys                         
import os
import numpy as np
import uncalled as unc

Opt = unc.config.Opt
MAPPER_OPTS = (
    Opt(("-t", "--threads"), ""),
    Opt("--num-channels", "read_buffer"),
    Opt(("-c", "--max-chunks"), "read_buffer"),
    Opt("--chunk-time", "read_buffer"),
    Opt("--config",
        type = str, 
        default = None, 
        required = False, 
        help = "Config file",
        dest = "_config_toml"
    ),
    Opt("--rna", fn="set_r94_rna")
)

def run(config):
    assert_exists(config.bwa_prefix + ".bwt")
    assert_exists(config.bwa_prefix + ".uncl")

    sys.stderr.flush()

    mapper = unc.MapPool(config)
    
    sys.stderr.write("Loading fast5s\n")
    for fast5 in load_fast5s(config.fast5_reader.fast5_files, config.fast5_reader.recursive):
        if fast5 != None:
            mapper.add_fast5(fast5)

    sys.stderr.flush()

    sys.stderr.write("Mapping\n")
    sys.stderr.flush()

    n = 0

    try:
        while mapper.running():
            t0 = time.time()
            for p in mapper.update():
                p.print_paf()
                n += 1
            dt = time.time() - t0;
            if dt < MAX_SLEEP:
                time.sleep(MAX_SLEEP - dt);
    except KeyboardInterrupt:
        pass
    
    sys.stderr.write("Finishing\n")
    mapper.stop()

OPTS = unc.index.BWA_OPTS + unc.fast5.FAST5_OPTS + MAPPER_OPTS

CMD = unc.config.Subcmd(
    "map", 
    "Map fast5 files to a DNA reference", 
    OPTS, run
)
