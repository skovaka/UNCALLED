#!/usr/bin/env python

import sys
import argparse
import numpy as np
from uncalled.sim_utils import SeqsumProfile
from uncalled.pafstats import parse_paf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates enrichment from read until simulation")
    parser.add_argument("-u", "--uncalled-fname", required=True, type=str, default=None, help="Simulator output PAF file")
    parser.add_argument("-c", "--cov-fname", required=True, type=str, default=None, help="BED file of control read coverage. Should be output from 'bedools intersect' of control read alignments and the target region(s)")
    parser.add_argument("-s", "--seq-sum", required=True, type=str, default=None, help="Control sequencing summary")
    parser.add_argument("-t", "--sim-speed", required=False, type=float, default=1, help="Speed that the simulator was run at")
    args = parser.parse_args()

    sys.stderr.write("Reading uncalled\n")
    sim_duration = 0
    unc_reads = dict()
    for p in parse_paf(args.uncalled_fname):

        end = (p.tags['st'][0]/4000) + (p.qr_len/450)
        sim_duration = max(sim_duration, end)

        v = (p.qr_len, p.tags.get('ej',(None,None))[0], p.tags.get('dl', (None,None))[0])

        if p.qr_name in unc_reads:
            unc_reads[p.qr_name].append(v)
        else:
            unc_reads[p.qr_name] = [v]

    sys.stderr.write("Reading minimap2\n")

    tgt_reads = dict()

    for line in open(args.cov_fname):
        ch,st,en,rd = line.split()[:4]
        tgt_reads[rd] = tgt_reads.get(rd,0)+int(en)-int(st)

    sys.stderr.write("Reading summary\n")

    read_info = dict()

    uc_counts = np.zeros((2,2))

    co = ct = 0
    uo = ut = 0

    ctl_reads = SeqsumProfile(args.seq_sum)
    ctl_reads.rm_scans()

    for i in range(len(ctl_reads)):

        read_id = ctl_reads.ids[i]
        start_time = ctl_reads.sts[i]
        duration = ctl_reads.lns[i]
        end_time = ctl_reads.ens[i]
        template_duration = ctl_reads.tds[i]
        seqlen = ctl_reads.bps[i]
        template_start = ctl_reads.tms[i]

        tgt_bp = tgt_reads.get(read_id, 0)

        end_time = start_time + duration

        ct += tgt_bp
        co += seqlen-tgt_bp

        unc_alns = unc_reads.get(read_id, None)

        if unc_alns == None:
            continue

        bpps = seqlen/template_duration

        for unc_est,eject_time,delay_time in unc_alns:
            ejected = eject_time != None

            if ejected:
                unclen = bpps * ((unc_est/450.0) + (delay_time/4000) - template_start + eject_time)
                if tgt_bp > 0: 
                    ut += unclen
                else:     
                    uo += unclen
            else:
                ut += tgt_bp
                uo += seqlen-tgt_bp

    up = 100*ut/(ut+uo)
    cp = 100*ct/(ct+co)

    co /= 1e6 
    ct /= 1e6 
    uo /= 1e6 
    ut /= 1e6 

    sys.stderr.write("\n")
    print("unc_on_bp\t%.6f" % (ut/args.sim_speed))
    print("unc_total_bp\t%.6f" % ((ut+uo)/args.sim_speed))
    print("cnt_on_bp\t%.6f" % ct)
    print("cnt_total_bp\t%.6f" % (ct+co))
