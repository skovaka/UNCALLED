#!/usr/bin/env python

import sys
import argparse
import numpy as np
from uncalled.sim_utils import SeqsumProfile
from uncalled.pafstats import parse_paf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates enrichment from read until simulation")
    parser.add_argument("-u", "--uncalled-fname", required=True, type=str, default=None, help="Simulator output PAF file")
    parser.add_argument("-s", "--seq-sum", required=True, type=str, default=None, help="Control sequencing summary")
    parser.add_argument("-m", "--minimap-fname", required=True, type=str, default=None, help="Minimap2 PAF file of the control reads aligned to a reference containing the target (or off-target, in the case of depletion) sequences")
    parser.add_argument("-x", "--bwa-prefix", required=True, type=str, default=None, help="BWA reference used during the simulation")
    parser.add_argument("--deplete", action='store_true')
    parser.add_argument("--enrich", action='store_true')
    parser.add_argument("-t", "--sim-speed", required=False, type=float, default=1, help="Speed that the simulator was run at")
    args = parser.parse_args()



    if args.deplete == args.enrich:
        sys.stderr.write("Need to specify one\n")
        sys.exit(1)

    sys.stderr.write("Reading reference annotation\n")
    ann_in = open("%s.ann" % args.bwa_prefix)
    nrefs = int(ann_in.readline().split()[1])
    ref_seqs = set()
    for i in range(nrefs):
        ref_seqs.add(ann_in.readline().split()[1])
        ann_in.readline()

    sys.stderr.write("Reading uncalled\n")
    sim_duration = 0
    unc_reads = dict()
    for p in parse_paf(args.uncalled_fname):

        end = (p.tags['st'][0]/4000) + (p.qr_len/450)
        sim_duration = max(sim_duration, end)

        v = (p.qr_len, p.tags.get('ej',(None,0))[0], p.tags.get('dl', (0,0))[0])
        if p.qr_name in unc_reads:
            unc_reads[p.qr_name].append(v)
        else:
            unc_reads[p.qr_name] = [v]

    sys.stderr.write("Reading minimap2\n")

    mm2_maps = [(p.qr_name, p.rf_name)
                for p in parse_paf(args.minimap_fname) 
                if p.is_mapped and p.tags['tp'][0] == 'P']

    mapped_reads = {q for q,r in mm2_maps}

    tgt_reads = {q for q,r in mm2_maps
                 if (args.deplete and r not in ref_seqs) or
                    (not args.deplete and r in ref_seqs)}

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

        ontgt = read_id in tgt_reads
        
        if ontgt: ct += seqlen
        else: co += seqlen

        unc_alns = unc_reads.get(read_id, None)

        if unc_alns == None:
            continue

        bpps = seqlen/template_duration

        for unc_est,eject_time,delay_time in unc_alns:
            ejected = eject_time != None

            if ejected:
                unclen = bpps * ((unc_est/450.0) + (delay_time/4000.0) + eject_time - template_start)
                if ontgt: 
                    ut += min(seqlen, unclen)
                else:     
                    uo += min(seqlen, unclen)

            elif ontgt:
                ut += seqlen
            else:         
                uo += seqlen

    co /= 1e6 
    ct /= 1e6 
    uo /= 1e6 
    ut /= 1e6 

    sys.stderr.write("\n")
    print("unc_on_bp\t%.6f" % (ut/args.sim_speed))
    print("unc_total_bp\t%.6f" % ((ut+uo)/args.sim_speed))
    print("cnt_on_bp\t%.6f" % ct)
    print("cnt_total_bp\t%.6f" % (ct+co))

