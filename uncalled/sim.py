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

import sys, os
import numpy as np
import argparse
from uncalled.pafstats import parse_paf
from time import time
from collections import Counter
import uncalled as unc

SAMP_RATE = 4000
CHS = np.arange(512)+1
NCHS = len(CHS)

NTICKS = 32
TICK_MOD = NCHS / NTICKS
PROG = " progress "
SP = '-'*( (NTICKS - len(PROG))//2 - 1)
PROG_HEADER = '|'+SP+PROG+SP+"|\n"

Opt = unc.ArgParser.Opt
OPTS = (
    Opt("fast5s", 
        nargs = '+', 
        type = str, 
        help = "Reads to unc. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
    ),
    Opt(("-r", "--recursive"), 
        action = "store_true"
    ),
    Opt("--ctl-seqsum", "simulator"),
    Opt("--unc-seqsum", "simulator"),
    Opt("--unc-paf",    "simulator"),
    Opt("--sim-speed",  "simulator"),
)

def run(conf):
    client = unc.Simulator(conf)
    load_sim(client, conf)

    for fast5 in load_fast5s(conf.fast5s, conf.recursive):
        if fast5 != None:
            client.add_fast5(fast5)

    client.load_fast5s()

    unc.realtime.run(conf, client)


def find_scans(sts,ens,mxs,max_block_gap=1,max_intv_gap=20,min_mux_frac=0.95):

    i = np.argsort(sts)
    sts = sts[i]
    ens = ens[i]
    mxs = mxs[i]

    blocks = list()
    bst = sts[0]
    ben = ens[0]

    for rst,ren in zip(sts[1:], ens[1:]):
        if rst - ben > max_block_gap:
            blocks.append( (bst,ben) )
            bst = rst
            ben = ren
        else:
            ben = max(ren,ben)
            #ben = ren
    blocks.append((bst,ben))

    scan_segs = list()
    scan = list()
    scan_gaps = list()


    prev_en = 0
    for bst,ben in blocks:

        if len(scan) > 0 and bst - scan[-1][1] > max_intv_gap:
            if len(scan) == 4:
                scan_segs.append(scan)
            scan = list()

        mux_counts = Counter(mxs[(sts >= bst) & (sts < ben)])
        mux_counts = [(c,m) for m,c in mux_counts.items()]
        top_count, top_mux = max(mux_counts)

        if top_count / sum((c for c,m in mux_counts)) >= min_mux_frac:

            if top_mux != 4 and len(scan) == 4:
                scan_segs.append(scan)

                scan_gaps.append((gap1,
                                  bst - scan[-1][1]))

                scan = list()

            if len(scan) > 0 and top_mux == len(scan):
                if ben - scan[-1][1] < max_intv_gap:
                    scan[-1] = (scan[-1][0],ben)

                elif top_mux == 1:
                    scan[0] = (bst,ben)
                    gap1 = bst - prev_en

            elif top_mux-1 == len(scan):
                scan.append( (bst,ben) )
                if len(scan) == 1:
                    gap1 = bst - prev_en

            else:
                scan = list()
        else:
            if len(scan) == 4:
                scan_segs.append(scan)
                scan_gaps.append((gap1,
                                  bst - scan[-1][1]))
            scan = list()

        prev_en = ben

    scans = list()

    for segs,gaps in zip(scan_segs, scan_gaps):
        scans.append((segs[0][0]-gaps[0], segs[-1][1]+gaps[1]))

    return scans

class SeqsumProfile:
    PROPS=['chs','sts','lns','mxs',
           'ids','ens','glns','gsts','tms','tds','bps']

    def __init__(self, fname, min_st=0, max_en=np.inf):
        infile = open(fname)
        header = infile.readline().split()
        ch_i = header.index("channel")
        st_i = header.index("start_time")
        ln_i = header.index("duration")
        mx_i = header.index("mux")
        id_i = header.index("read_id")
        tm_i = header.index("template_start")
        td_i = header.index("template_duration")
        bp_i = header.index("sequence_length_template")

        ids=list()
        chs=list()
        sts=list()
        mxs=list() 
        lns=list()
        ens=list()
        tms=list()
        tds=list()
        bps=list()

        SIZE = os.path.getsize(fname)
        MOD = SIZE / NTICKS
        bts = 0
        sys.stderr.write("=")

        for line in infile:
            bts += len(line)
            if bts > MOD:
                bts = 0
                sys.stderr.write("=")
                sys.stderr.flush()

            tabs = line.split()
            st = float(tabs[st_i])
            ln = float(tabs[ln_i])
            en = st+ln
            if st < min_st or en > max_en: continue

            sts.append(st)
            lns.append(ln)
            ens.append(en)
            chs.append(int(tabs[ch_i]))
            mxs.append(int(tabs[mx_i]))
            ids.append(tabs[id_i])
            tms.append(float(tabs[tm_i])-st)
            tds.append(float(tabs[td_i]))
            bps.append(int(tabs[bp_i]))

        sys.stderr.write("\n")
        sys.stderr.write("Procesing run...................\n")

        ids,chs,sts,mxs,lns,ens,tms,tds,bps=map(np.array, (ids,chs,sts,mxs,lns,ens,tms,tds,bps))

        self.ids,self.chs,self.sts,self.mxs,self.lns,self.ens,self.tms,self.tds,self.bps = ids,chs,sts,mxs,lns,ens,tms,tds,bps

        self.sort(np.argsort(sts))

        self.chodr = CHS
        self.chcts = np.array([np.sum(self.chs==ch) 
                               for ch in CHS])

        self.duraiton = np.max(ens)

    def rm_scans(self):
        scans = find_scans(self.sts, self.ens, self.mxs)
        bounds = list()

        sh = 0
        for st,en in (scans):
            m = np.flatnonzero((self.sts+sh >= st) & (self.ens+sh <= en))

            for pr in SeqsumProfile.PROPS:
                a = getattr(self, pr, None)
                if a is not None:
                    setattr(self, pr, np.delete(a,m))

            l = en-st
            bounds.append(st-sh)

            self.sts[self.sts+sh >= st] -= l
            self.ens[self.ens+sh >= st] -= l
            sh += l

        bounds.append(np.max(self.ens))

        self.chcts = np.array([np.sum(self.chs==ch) 
                               for ch in CHS])

        return np.array(bounds)
    
    def compute_eject_delays(self, fname):
        self.dls = np.zeros(len(self.sts))
        self.dls.fill(np.inf)

        idxs = {self.ids[i] : i for i in range(len(self.ids))}
        tlns = self.lns - self.tms

        for p in parse_paf(open(fname)):
            i = idxs.get(p.qr_name, None)
            if i != None and ('ej' in p.tags or 'ub' in p.tags):
                ej = (p.tags['ej'] if 'ej' in p.tags else p.tags['ub'])[0]
                self.dls[i] = max(0, tlns[i] - ((p.qr_len/450.0)+ej))

    def compute_gaps(self):
        self.gsts = np.zeros(len(self.ids))
        self.glns = np.zeros(len(self.ids))

        for ch in CHS:
            cids = self.ids[self.chs==ch]
            csts = self.sts[self.chs==ch]
            cens = self.ens[self.chs==ch]

            gsts = np.insert(cens[:-1], 0, 0)
            glns = csts - gsts

            self.gsts[self.chs==ch] = gsts
            self.glns[self.chs==ch] = glns

    def chsort(self, odr):
        self.chodr = self.chodr[odr]
        self.chcts = self.chcts[odr]


    def sort(self, order=None):
        for pr in SeqsumProfile.PROPS:
            a = getattr(self, pr, None)
            if a is not None:
                setattr(self, pr, a[order])

    def __len__(self):
        return len(self.sts)

def sec_to_samp(sec, coef=1.0):
    return int(np.round(sec*SAMP_RATE*coef))

#def write_itv(out, ch, sc, st_sec, en_sec, time_scale):
#    st_samp = np.round(st_sec*SAMP_RATE*time_scale)
#    en_samp = np.round(en_sec*SAMP_RATE*time_scale)
#    itvs_out.write("%d\t%d\t%d\t%d\n" % (ch,sc,st_samp,en_samp))
#
#def write_gap(out, ch, sc, ln_sec):
#    ln_samp = np.round(ln_sec*SAMP_RATE)
#    out.write("%d\t%d\t%d\n" % (ch, sc, ln_samp))

def load_sim(client, conf):

    t0 = time()

    sys.stderr.write("Loading UNCALLED PAF............\n")
    unc = SeqsumProfile(conf.unc_seqsum)

    unc_scans = unc.rm_scans()
    unc.compute_gaps()

    unc.compute_eject_delays(conf.unc_paf)
    delays = unc.dls[unc.dls != np.inf]
    DELAY = np.median(delays)
    unc.chsort(np.argsort(unc.chcts))

    sys.stderr.write("Generating pattern..............\n")

    ACTIVE_THRESH = np.median(unc.glns)+np.std(unc.glns)

    for ch in CHS:

        if (ch-1) % TICK_MOD == 0:
            sys.stderr.write("=")
            sys.stderr.flush()

        ch_i = unc.chs==ch

        if not np.any(ch_i):
            continue

        gsts = unc.gsts[ch_i]
        glns = unc.glns[ch_i]

        ch_intervals = list()

        #break scan intervals at long gaps
        sc = 0
        itv_st = 0
        for br in np.flatnonzero(glns >= ACTIVE_THRESH):
            act_en = gsts[br] #find first long gap (break)
            #print(ch, itv_st, _en)

            #add all full intervals preceding break
            while unc_scans[sc+1] < act_en:
                itv_en = conf.scan_intv_time

                st_samp = sec_to_samp(itv_st-unc_scans[sc], conf.sim_speed)
                en_samp = sec_to_samp(itv_en, conf.sim_speed)
                client.add_intv(ch,sc,st_samp,en_samp)

                #write_itv(itvs_out, ch, sc, itv_st-unc_scans[sc], itv_en, conf.sim_speed)

                itv_st = unc_scans[sc+1]
                sc += 1

            #add partial intervals before break
            if itv_st != act_en:
                #if np.any((glns < ACTIVE_THRESH) & (gsts > itv_st) & (gsts < act_en)):
                st_samp = sec_to_samp(itv_st-unc_scans[sc], conf.sim_speed)
                en_samp = sec_to_samp(act_en-unc_scans[sc], conf.sim_speed)
                client.add_intv(ch,sc,st_samp,en_samp)

                #write_itv(itvs_out, ch, sc, itv_st-unc_scans[sc], act_en-unc_scans[sc], conf.sim_speed)

            itv_st = act_en + glns[br]
            
            #skip intervals before interval starts
            while unc_scans[sc+1] < itv_st:
                sc += 1

        last = np.max(unc.ens[ch_i]) #time of last read

        #add intervals between last break and final read
        while sc < len(unc_scans)-1 and unc_scans[sc] < last:
            itv_en = min(last - unc_scans[sc], conf.scan_intv_time)

            st_samp = sec_to_samp(itv_st-unc_scans[sc], conf.sim_speed)
            en_samp = sec_to_samp(itv_en, conf.sim_speed)

            client.add_intv(ch,sc,st_samp,en_samp)
            #write_itv(itvs_out,ch,sc,itv_st-unc_scans[sc],itv_en, conf.sim_speed)

            itv_st = unc_scans[sc+1]
            sc += 1

        #write short gaps within each interval
        for sc in range(len(unc_scans)-1):
            sc_i = (gsts > unc_scans[sc]) & ((gsts+glns) <= unc_scans[sc+1])
            for ln in glns[sc_i]:
                if ln < ACTIVE_THRESH and ln > 0:
                    client.add_gap(ch, sc, sec_to_samp(ln))
                    #write_gap(gaps_out, ch, sc, ln)

            for dl in unc.dls[ch_i][sc_i]:
                if dl != np.inf:
                    #write_gap(delays_out, ch, sc, dl)
                    #write_gap(delays_out, ch, sc, DELAY)
                    client.add_delay(ch, sc, sec_to_samp(DELAY))


    sys.stderr.write("\n\n")

    #gaps_out.close()
    #itvs_out.close()
    #delays_out.close()

    sys.stderr.write("Loading control PAF.............\n")
    ctl = SeqsumProfile(conf.ctl_seqsum)
    ctl.rm_scans()
    ctl.chsort(np.argsort(ctl.chcts))

    sys.stderr.write("Ordering reads..................\n")
    #sys.stderr.write(PROG_HEADER)

    #Channels with any reads recieve minimum read count
    min_const = np.zeros(NCHS)
    min_const[unc.chcts > 0] = conf.min_ch_reads

    tgt_total = np.sum(ctl.chcts)

    #Max reads proportional to unc read counts
    max_prpl = tgt_total * unc.chcts / np.sum(unc.chcts)

    #Combine and min const with re-scaled proportional to fit max reads
    remain = max_prpl - min_const
    remain_clp = np.clip(remain, 0, np.inf)
    tgt_counts = min_const + (np.sum(remain) * remain_clp / np.sum(remain_clp))

    #Round and adjust for error
    tgt_counts = np.round(tgt_counts).astype(int)
    dr = -1 if np.sum(tgt_counts) > tgt_total else 1
    i = len(tgt_counts)-1
    while np.sum(tgt_counts) != tgt_total:
        tgt_counts[i] += dr
        i -= 1

    diff = ctl.chcts - tgt_counts

    odr = np.flip(np.argsort(diff),0)
    diff = diff[odr]
    tgt_counts = tgt_counts[odr]
    ctl.chsort(odr)
    unc.chsort(odr)

    sim_reads = [None for c in CHS]
    extra = list()
    e = 0

    for i in range(NCHS):
        if i % TICK_MOD == 0:
            sys.stderr.write("=")
            sys.stderr.flush()
        
        j = ctl.chs==ctl.chodr[i]
        
        ctl_reads = list(zip(ctl.ids[j], ctl.tms[j]))
        new_reads = None

        tgt = tgt_counts[i]

        if diff[i] >= 0:
            new_reads = ctl_reads[:tgt]
            
            if diff[i] > 0:
                extra.append(ctl_reads[tgt:])
        else:
            if e >= len(extra):
                sys.stderr.write("Not enough reads? maybe should haved checked earlier\n")
                sys.exit(1)

            new_reads = ctl_reads
            while len(new_reads) < tgt and e < len(extra):
                needed = tgt-len(new_reads)
                if len(extra[e]) > needed:
                    new_reads += extra[e][:needed]
                    extra[e] = extra[e][needed:]
                else:
                    new_reads += extra[e]
                    e += 1

            if len(new_reads) < tgt:
                sys.stderr.write("Not enough reads? again? not sure I should be here\n")
                sys.exit(1)

        sim_reads[unc.chodr[i]-1] = new_reads

    sys.stderr.write("\n")

    for ch in CHS:
        for rd,tm in sim_reads[ch-1]:
            client.add_read(ch, rd, sec_to_samp(tm))
            #reads_out.write("%3d %s %d\n" % (ch, rd, np.round(tm*SAMP_RATE)))
            #reads_out.write("%3d %s %d\n" % (ch, rd, 0))
    #reads_out.close()
