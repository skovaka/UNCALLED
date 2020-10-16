#!/usr/bin/env python

import sys
import numpy as np
import re
import argparse

class PafEntry:
    def __init__(self, line, tags=None):
        fromstr = type(line) != list
        if fromstr:
            tabs = line.split()
        else:
            tabs = line

        self.qr_name = tabs[0]
        self.qr_len = int(tabs[1])
        self.is_mapped = tabs[4] != ("*" if fromstr else None)

        if self.is_mapped:
            self.qr_st = int(tabs[2])
            self.qr_en = int(tabs[3])
            self.is_fwd = tabs[4] == ('+' if fromstr else True)
            self.rf_name = tabs[5]
            self.rf_len = int(tabs[6])
            self.rf_st = int(tabs[7])
            self.rf_en = int(tabs[8])
            self.match_num = int(tabs[9])
            self.aln_len = int(tabs[10])
            self.qual = int(tabs[11])
        else:
            self.qr_st = 1
            self.qr_en = self.qr_len
            self.is_fwd=self.rf_name=self.rf_len=self.rf_st=self.rf_en=self.match_num=self.aln_len=self.qual=None

        self.tags = dict() if tags==None else tags 
        for k,t,v in (s.split(":") for s in tabs[12:]):
            if t == 'f':
                v = float(v)
            elif t == 'i':
                v = int(v)
            elif t not in ['A', 'Z','B', 'H']:
                sys.stderr.write("Error: invalid tag type \"%s\"\n" % t)
                sys.exit(1)
            self.tags[k] = (v,t)

    def rev(self):
        return PafEntry( [self.rf_name, self.rf_len, self.rf_st, self.rf_en, self.is_fwd, 
                          self.qr_name, self.qr_len, self.qr_st, self.qr_en, self.match_num,
                          self.aln_len, self.qual], self.tags )

    def get_tag(self, k):
        return self.tags.get(k, (None,None))[0]
    
    def set_tag(self, k, v, t=None):
        if t == None:
            if isinstance(v, int):
                t = 'i'
            elif isinstance(v, float):
                t = 'f'
            else:
                t = 'Z'
        self.tags[k] = (v,t)

    def qry_loc(self):
        return (self.qr_name, self.qr_st, self.qr_en)

    def ref_loc(self):
        return (self.rf_name, self.rf_st, self.rf_en)

    def ext_ref(self, ext=1.0):
        st_shift = int(self.qr_st*ext)
        en_shift = int((self.qr_len - self.qr_en)*ext)

        if self.is_fwd:
            return (max(1, self.rf_st - st_shift),
                    min(self.rf_len, self.rf_en + en_shift))
        else:
            return (max(1, self.rf_st - en_shift),
                    min(self.rf_len, self.rf_en + st_shift))
    
    def overlaps(self, paf2, ext=0.0):
        st1, en1 = self.ext_ref(ext)
        st2, en2 = paf2.ext_ref(ext)
        return (self.is_mapped and paf2.is_mapped and self.rf_name.startswith(paf2.rf_name) and 
                max(st1, st2) <= min(en1, en2))

    def contains(self, paf2):
        return (self.is_mapped and paf2.is_mapped and self.rf_name.startswith(paf2.rf_name) and 
                self.rf_st <= paf2.rf_st and self.rf_en >= paf2.rf_en)

    def __lt__(self, paf2):
        return self.qr_name < self.qr_name

    def __str__(self):
        tagstr = "\t".join( (":".join([k,v[1],str(v[0])]) for k,v in self.tags.items()))
        if self.is_mapped:
            s = "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s" % (
                 self.qr_name, self.qr_len, self.qr_st, self.qr_en, 
                 '+' if self.is_fwd else '-', self.rf_name, 
                 self.rf_len, self.rf_st, self.rf_en, 
                 self.match_num, self.aln_len, self.qual, tagstr)
        else:
            s = "\t".join((self.qr_name,str(self.qr_len)) + ("*",)*10 + (tagstr,))

        return s


def parse_paf(infile, max_load=None):
    if isinstance(infile, str):
        infile = open(infile)
    c = 0
    for l in infile:
        if l[0] == "#": continue
        if max_load != None and c >= max_load: break
        yield PafEntry(l)
        c += 1

def paf_ref_compare(qry, ref, ret_qry=True, check_locs=True, ext=1.5):
    if type(ref) == dict:
        ref_locs = ref
    else:
        ref_locs = dict()
        for r in ref:
            l = ref_locs.get(r.qr_name, None)
            if l == None:
                ref_locs[r.qr_name] = [r]
            else:
                l.append(r)

    tp = list()
    tn = list()
    fp = list()
    fn = list()
    fp_unmap = list()

    for q in qry:
        rs = ref_locs.get(q.qr_name, [None])
        if q.is_mapped:
            if rs == [None] or not rs[0].is_mapped:
                fp_unmap.append(q if ret_qry else rs[0])
                continue
            match = False
            for r in rs:
                if ((check_locs and q.overlaps(r, ext)) or 
                    (not check_locs and q.rf_name == r.rf_name)):
                        match = True
                        tp.append(q if ret_qry else r)
                        break
            if not match:
                fp.append(q if ret_qry else r)
        else:
            if rs == [None] or not rs[0].is_mapped:
                tn.append(q if ret_qry else rs[0])
            else:
                fn.append(q if ret_qry else rs[0])

    return tp, tn, fp, fn, fp_unmap

def add_opts(parser):
    parser.add_argument("infile", type=str, help="PAF file output by UNCALLED")
    parser.add_argument("-n", "--max-reads", required=False, type=int, default=None, help="Will only look at first n reads if specified")
    parser.add_argument("-r", "--ref-paf", required=False, type=str, default=None, help="Reference PAF file. Will output percent true/false positives/negatives with respect to reference. Reads not mapped in reference PAF will be classified as NA.")
    parser.add_argument("-a", "--annotate", action='store_true', help="Should be used with --ref-paf. Will output an annotated version of the input with T/P F/P specified in an 'rf' tag")

def run(args):
    locs = [p for p in parse_paf(args.infile, args.max_reads)]

    num_mapped = sum([p.is_mapped for p in locs])

    statsout = sys.stderr if args.annotate else sys.stdout

    statsout.write("Summary: %d reads, %d mapped (%.2f%%)\n\n" % (len(locs), num_mapped, 100*num_mapped/len(locs)))

    if args.ref_paf != None:
        statsout.write("Comparing to reference PAF\n")
        tp, tn, fp, fn, fp_unmap = paf_ref_compare(locs, parse_paf(args.ref_paf))
        ntp,ntn,nfp,nfn,nfp_unmap = map(len, [tp, tn, fp, fn, fp_unmap])
        n = len(locs)

        statsout.write("     P     N\n")
        statsout.write("T %6.2f %5.2f\n" % (100*ntp/n, 100*ntn/n))
        statsout.write("F %6.2f %5.2f\n" % (100*(nfp)/n, 100*nfn/n))
        statsout.write("NA: %.2f\n\n" % (100*nfp_unmap/n))

        if args.annotate:
            group_labels = [(tp, "tp"), 
                            (tn, "tn"), 
                            (fp, "fp"),
                            (fn, "fn"),
                            (fp_unmap, "na")]

            for grp,lab in group_labels:
                for p in grp:
                    p.set_tag("rf", lab, "Z")
                    sys.stdout.write("%s\n" % p)

    if locs[0].get_tag('mt') != None:
        map_ms = np.array([p.get_tag('mt') for p in locs if p.is_mapped])
        map_bp = np.array([p.qr_en for p in locs if p.is_mapped])
        map_bpps = 1000*map_bp/map_ms

        statsout.write("Speed            Mean    Median\n")
        statsout.write("BP per sec: %9.2f %9.2f\n" % (np.mean(map_bpps), np.median(map_bpps)))
        statsout.write("BP mapped:  %9.2f %9.2f\n" % (np.mean(map_bp),   np.median(map_bp)))
        statsout.write("MS to map:  %9.2f %9.2f\n" % (np.mean(map_ms),   np.median(map_ms)))
