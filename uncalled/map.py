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

Opt = unc.ArgParser.Opt
OPTS = (
    Opt(("-t", "--threads"), ""),
    Opt("--num-channels", "read_buffer"),
    Opt(("-e", "--max-events"), "mapper"),
    Opt(("-c", "--max-chunks"), "read_buffer"),
    Opt("--chunk-time", "read_buffer"),
    Opt("--conf", 
        type = str, 
        default = None, 
        required = False, 
        help = "Config file"
    ),
    Opt("--rna", 
        action = "store_true",
        help = "Will use RNA parameters if set"
    ),

    #TODO move to different parser set
    Opt(("-o", "--out-prefix"), 
        type = str, 
        default = None, 
        required = False, 
        help = "Output prefix"
    ),
    Opt("--mm2", 
        type = str, 
        default = None, 
        required = False, 
        help = "Minimap2 PAF file for comparison"
    ),
) #end MAP_OPTS

def run(conf):
    assert_exists(conf.bwa_prefix + ".bwt")
    assert_exists(conf.bwa_prefix + ".uncl")

    sys.stderr.flush()

    mapper = unc.MapPool(conf)
    
    sys.stderr.write("Loading fast5s\n")
    for fast5 in load_fast5s(conf.fast5_reader.fast5_files, conf.fast5_reader.recursive):
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
