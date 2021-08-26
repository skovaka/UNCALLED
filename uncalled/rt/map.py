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

import sys            
import time

from ..index import check_prefix
from ..argparse import Opt, CONFIG_PARAM
from ..index import BWA_OPTS
from ..fast5 import FAST5_OPTS
from _uncalled import MapPoolOrd

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
        dest = CONFIG_PARAM
    ),
    Opt("--rna", fn="set_r94_rna")
)

OPTS = BWA_OPTS + FAST5_OPTS + MAPPER_OPTS

def main(conf):
    """Map fast5 files to a DNA reference"""
    check_prefix(conf.bwa_prefix)

    mapper = MapPoolOrd(conf)
    
    sys.stderr.write("Loading fast5s\n")
    sys.stderr.flush()
    mapper.load_fast5s()

    sys.stderr.write("Mapping\n")
    sys.stderr.flush()

    MAX_SLEEP = 0.01
    try:
        while mapper.running():
            t0 = time.time()
            for ch,nm,paf in mapper.update():
                paf.print_paf()

            dt = time.time() - t0;
            if dt < MAX_SLEEP:
                time.sleep(MAX_SLEEP - dt);
    except KeyboardInterrupt:
        pass
    
    sys.stderr.write("Finishing\n")
    mapper.stop()

