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

from .. import index
from .map import MAPPER_OPTS
from ..argparse import Opt, MutexOpts
from .read_until import ReadUntilClient


REALTIME_OPTS = (
    MutexOpts("realtime_mode", [
        Opt(("-D", "--deplete"), fn="set_rt_deplete"),
        Opt(("-E", "--enrich"), fn="set_rt_enrich"),
    ]),
    MutexOpts("active_chs", [
        Opt("--even", fn="set_active_chs_even", help="world"),
        Opt("--odd", fn="set_active_chs_odd", help="Hello"),
    ])
)

OPTS = index.BWA_OPTS + MAPPER_OPTS + (
    Opt("--host", "realtime"),
    Opt("--port", "realtime"),
    Opt("--duration", "realtime"),
) + REALTIME_OPTS

def main(config, client=None):
    """Perform real-time targeted sequencing"""

    index.check_prefix(config.bwa_prefix)

    pool = None
    client = None

    #TODO make simulator indisinguishable from client, maybe w/ python wrapper
    if client is None:
        sim = False
        client = ReadUntilClient(config.host, config.port, config.chunk_time, config.num_channels)
    else:
        sim = True

    try:
        if not client.run():
            sys.stderr.write("Error: failed to run client\n")
            sys.exit(1)

        deplete = config.realtime_mode == unc.RealtimePool.DEPLETE
        even = config.active_chs == unc.RealtimePool.EVEN #TODO: do within mapper

        if not sim:
            raw_type = str(client.signal_dtype)

        pool = unc.RealtimePool(config)

        chunk_times = [time.time() for c in range(config.num_channels)]
        unblocked = [None for c in range(config.num_channels)]

        if config.duration == None or config.duration == 0:
            end_time = float("inf")
        else:
            end_time = config.duration*60*60

        while client.is_running:
            t0 = time.time()

            for ch, nm, paf in pool.update():
                t = time.time()-chunk_times[ch-1]
                if paf.is_ended():
                    paf.set_float(unc.Paf.ENDED, t)
                    client.stop_receiving_read(ch, nm)

                elif (paf.is_mapped() and deplete) or not (paf.is_mapped() or deplete):

                    if sim or client.should_eject():
                        paf.set_float(unc.Paf.EJECT, t)
                        u = client.unblock_read(ch, nm)

                        if sim:
                            paf.set_int(unc.Paf.DELAY, u)

                        unblocked[ch-1] = nm
                    else:
                        paf.set_float(unc.Paf.IN_SCAN, t)
                        client.stop_receiving_read(ch, nm)

                else:
                    paf.set_float(unc.Paf.KEEP, t)
                    client.stop_receiving_read(ch, nm)

                paf.print_paf()

            if sim:
                read_batch = client.get_read_chunks()
                for channel, read in read_batch:
                    if even and channel % 2 == 1:
                        client.stop_receiving_read(channel, read.number)
                    else:
                        if unblocked[channel-1] == read.number:
                            sys.stdout.write("# recieved chunk from %s after unblocking\n" % read.id)
                            continue

                        chunk_times[channel-1] = time.time()
                        pool.add_chunk(read)
       
            else:

                read_batch = client.get_read_chunks(batch_size=client.queue_length)
                for channel, read in read_batch:
                    if even and channel % 2 == 1:
                        client.stop_receiving_read(channel, read.number)
                    else:
                        if unblocked[channel-1] == read.number:
                            sys.stdout.write("# recieved chunk from %s after unblocking\n" % read.id)
                            continue

                        chunk_times[channel-1] = time.time()
                        pool.add_chunk(unc.Chunk(read.id, 
                                                     channel, 
                                                     read.number,
                                                     read.chunk_start_sample,
                                                     raw_type,
                                                     read.raw_data))


            if client.get_runtime() >= end_time:
                if not sim:
                    client.reset()
                client = None
                break

            dt = time.time() - t0;
            if dt < MAX_SLEEP:
                time.sleep(MAX_SLEEP - dt);

    except KeyboardInterrupt:
        sys.stderr.write("Keyboard interrupt\n")

    except Exception as e:
        sys.stderr.write(traceback.format_exc())

    #client.log("Finished")

    if client != None and not sim:
        client.reset()

    if pool != None:
        pool.stop_all()


