#!/usr/bin/env python

import sys, os
import numpy as np
import argparse
import bisect
from collections import defaultdict
import re

SEED_LEN = 22
MIN_CLUST = 25
CHUNK_LEN = 4000

IS_RNA = False#True

if IS_RNA:
    SAMPLE_RATE = 3012
    BCE_STRIDE = 10
else:
    SAMPLE_RATE = 4000
    BCE_STRIDE = 5

#Guppy basecalled event stride and kmer length
BCE_K = 4

MAX_CHUNK_DEF = 3

#Cigar parsing
CIG_OPS_STR = "MIDNSHP=X"
CIG_RE = re.compile("(\d+)(["+CIG_OPS_STR+"])")
CIG_OPS = set(CIG_OPS_STR)
CIG_INCR_ALL = {'M','=', 'X'}
CIG_INCR_RD = CIG_INCR_ALL | {'I','S'}
CIG_INCR_RF = CIG_INCR_ALL | {'D','N'}

BASES = "ACGT"

def coords_inverted(fwd):
    return not (fwd or IS_RNA) or (fwd and IS_RNA)

class DebugParser:
    CONF_PAD_COEF = 2

    def __init__(self, 
                 debug_prefix, 
                 unc_paf, 
                 mm2_paf=None,
                 bce_moves=None,
                 bwa_index=None,
                 min_chunk=None,
                 max_chunk=None,
                 min_samp=None,
                 max_samp=None,
                 tgt_cid=None,
                 max_path_fm=0,
                 load_seeds=True,
                 load_events=True,
                 load_paths=True,
                 load_conf=True):

        self.rid = unc_paf.qr_name

        #Cofident seed cluster ID, confident event, and reference length
        #(conf_evt is where mapping would normally end)
        self.conf_cid = unc_paf.tags.get('sc', (None,)*2)[0]
        self.conf_evt = unc_paf.tags.get('ce', (None,)*2)[0]
        self.conf_samp = None

        self.mm2_paf = mm2_paf

        self.idx = bwa_index

        #if self.tgt_cid is None:
        #    self.tgt_cid = self.conf_cid

        if min_samp is None and min_chunk is None:
            self.min_samp = self.min_chunk = 0
        elif min_chunk is None:
            self.min_samp = min_samp
            self.min_chunk = (min_samp-1) // CHUNK_LEN
        elif min_samp is None:
            self.min_samp = self.min_chunk * SAMPLE_RATE
            self.min_chunk = min_chunk
        else:
            self.min_samp = max(min_samp, min_chunk * SAMPLE_RATE)
            self.min_chunk = max(min_chunk, (min_samp-1) // CHUNK_LEN)


        if max_samp is None and max_chunk is None:
            if self.conf_evt is None:
                self.max_chunk = MAX_CHUNK_DEF
                self.max_samp = self.max_chunk*SAMPLE_RATE
            else:
                self.max_chunk = self.max_samp = None
        elif max_samp is None:
            self.max_chunk = max_chunk
            self.max_samp = max_chunk * SAMPLE_RATE
        elif max_chunk is None:
            self.max_chunk = max(0, (max_samp-1) // CHUNK_LEN)
            self.max_samp = max_samp
        else:
            self.max_chunk = max(max_chunk, (max_samp-1) // CHUNK_LEN)
            self.max_samp = max(max_samp, self.max_chunk * SAMPLE_RATE)

        self.min_evt = 0 if self.min_samp == 0 else None
        self.max_evt = None

        self.chunk_evt_bounds = [0]

        #TODO below should be w/r/t reference coordinates?
        #allow argument specified
        #or detect based on confident cluster (or just specify seed cluster)

        #TBD in parse_seeds
        self.conf_len = None 
        self.conf_clust = None
        self.max_clust = None

        self.max_clust_len = None

        if mm2_paf is not None:
            self.ref_name = mm2_paf.rf_name
            self.fwd = mm2_paf.is_fwd

            self.inv_coords = coords_inverted(mm2_paf.is_fwd)
            if self.inv_coords:
                self.min_ref = -mm2_paf.rf_en
            else:
                self.min_ref = mm2_paf.rf_st

        else:
            self.inv_coords = None
            self.min_ref = None
            self.ref_name = None

        self.max_ref = None

        #TODO clean these up
        self.conf_pbs = dict()
        self.conf_dots = set()
        self.path_lens = defaultdict(int)
        self.seed_kmers = dict()
        self.conf_ref_bounds = None

        self.evts_loaded = False
        if load_events:
            self.evts_in = open(debug_prefix + self.rid + "_events.tsv")
            self.parse_events()

        self.seeds_loaded = False
        if load_seeds:
            self.seeds_in  = open(debug_prefix + self.rid + "_seeds.bed")
            self.parse_seeds()

        self.bc_loaded = False
        if mm2_paf is not None and bce_moves is not None:
            self.parse_bc_aln(bce_moves)
        else:
            print(mm2_paf, bce_moves)

        #TODO: could look at path bufs, but need to specify coords
        #also can store max cluster, but need to store cids
        self.empty = self.mm2_paf is None and self.conf_clust is None
        if self.empty:
            return

        self.max_path_fm = max_path_fm
        if max_path_fm > 0:
            if self.idx is None:
                sys.stderr.write("Error: must include BWA index to include paths with FM overlaps\n")
                sys.exit(1)

            #fwd =  (mm2_paf is not None and mm2_paf.is_fwd or
            #        mm2_paf is None and self.max_clust.fwd)

            #inv = coords_inverted(fwd)
            if not self.inv_coords:
                st = self.min_ref
                en = self.max_ref
            else:
                st = -self.max_ref
                en = -self.min_ref

            ############
            #sa_st = self.idx.get_sa_loc(self.ref_name, st)

            #seq = ""
            #for i in range(en-st+1):
            #    c = BASES[self.idx.get_base(sa_st+i)]
            #    print(i, c)
            #    seq = seq + c

            #print(seq)
            ############

            fwd_fms, rev_fms = self.idx.range_to_fms(
                    self.ref_name, st, en+1
            )
            #fms = fwd_fms if fwd else rev_fms
            fms = fwd_fms if not self.inv_coords else rev_fms

            #print(fms)

            self.fm_to_ref = {fms[i] : i for i in range(len(fms))}

            fm_ord = np.argsort(fms)
            if self.bc_loaded:
                #print(len(fms), self.bce_refs[-1]-self.bce_refs[0], self.bce_refs[0])

                evt = np.searchsorted(self.evt_ens, self.bce_samps[0])
                evt_range = range(0, len(self.events))
                ref_range = range(self.bce_refs[0], self.bce_refs[-1]+1)

                self.range_fm_evts = [set() for r in range(len(fms))]

                for bce_s,bce_r in zip(self.bce_samps, self.bce_refs):
                    while evt in evt_range and self.evt_ens[evt] < bce_s:
                        evt += 1
                    if evt > len(self.evt_ens): break

                    N = 1
                    for e in range(evt-N, evt+N+1):
                        for r in range(bce_r-N, bce_r+N+1):
                            if e in evt_range and r in ref_range:
                                self.range_fm_evts[r+3].add(e)

                self.range_fm_evts = np.array(self.range_fm_evts)[fm_ord]

            self.range_fms = np.array(fms)[fm_ord]

        else:
            self.fwd_fms = self.rev_fms = None

        #TODO: paths
        self.paths_loaded = False
        if load_paths:
            self.paths_fname = debug_prefix + self.rid + "_paths.tsv"
            self.parse_paths()

        self.conf_loaded = False
        if load_conf:
            self.conf_fname = debug_prefix + self.rid + "_conf.tsv"
            self.parse_conf()

    def parse_events(self, incl_norm=True):
        if self.evts_loaded: return False
        self.evts_in.readline()

        self.events = list()
        self.norms = list()

        self.win_means = list()
        self.win_stdvs = list()
        self.evt_mask = list()

        self.evt_mask_map = list()

        unmask_evts = 0

        evt = 0
        chunk = 0
        next_chunk_samp = CHUNK_LEN

        for line in self.evts_in:
            tabs = line.split()
            st,ln,mask = map(int, tabs[:2] + tabs[-1:])
            mn,sd,norm_sc,norm_sh,win_mn,win_sd = map(float, tabs[2:-1])
            en = st+ln

            if st >= next_chunk_samp:
                chunk += 1
                next_chunk_samp += CHUNK_LEN
                self.chunk_evt_bounds.append(evt)

            if evt == self.conf_evt:
                self.conf_samp = st

                if self.max_samp is None:
                    self.max_samp = (chunk+1)*CHUNK_LEN-1
                    self.max_chunk = chunk

            if en < self.min_samp: 
                evt += 1
                continue

            if self.min_evt is None:
                self.min_evt = evt

            if self.max_samp is not None and st > self.max_samp: 
                break
            
            self.events.append( (st,ln,mn,sd) )
            self.norms.append( (norm_sc,norm_sh) )

            self.win_means.append(win_mn)
            self.win_stdvs.append(win_sd)
            self.evt_mask.append(mask == 1)

            if self.evt_mask[-1]:
                self.evt_mask_map.append(evt)

            evt += 1

        if len(self.events) == 0:
            self.evts_loaded = False
            return False

        if self.max_evt is None or self.max_evt >= evt: 
            self.max_evt = evt

        if self.max_samp is None or st+ln < self.max_samp:
            self.max_samp = st+ln
            self.max_chunk = (self.max_samp-1) // CHUNK_LEN

        self.evts_loaded = True

        self.evt_sts,self.evt_lns,self.evt_mns,self.evt_sds = map(np.array, zip(*self.events))
        self.evt_ens = self.evt_sts + self.evt_lns

        self.win_means = np.array(self.win_means)
        self.win_stdvs = np.array(self.win_stdvs)
        self.evt_mask = np.array(self.evt_mask)
        self.evt_mask_map = np.array(self.evt_mask_map)

        return True
    
    def normed_event(self, e):
        scale,shift = self.norms[e]
        return scale*self.events[e][2]+shift

    def parse_seeds(self, expire_coef=None):
        if self.seeds_loaded: return False

        if expire_coef != None:
            SeedCluster.EXPIRE_COEF = expire_coef

        self.clusts = dict()
        self.clusts_exp = list()
    
        clust_ids = set()

        conf_clust = None
        for line in self.seeds_in:
            rf,st,en,name,strand = line.split()
            st,en,evt,pb,cid = map(int, [st,en]+name.split(":"))

            #handle expiration
            #should this (still) be here?
            clust = self.clusts.get(cid, None)
            replace = (
                clust is not None and 
                clust.expired(evt) and
                (cid != self.conf_cid or evt < self.conf_evt)
            )

            if clust == None or replace:
                if replace:
                    self.clusts_exp.append(clust)

                clust = SeedCluster(evt,rf,st,en,strand=="+",cid)
                self.clusts[cid] = clust

                clust_ids.add(cid)

            else:
                clust.add_seed(evt,rf,st,en)

            if clust.id == self.conf_cid:
                conf_clust = clust
                if not coords_inverted(clust.fwd):
                    #self.conf_pbs[(evt,pb)] = en
                    self.conf_pbs[(evt,pb)] = st
                else:
                    #self.conf_pbs[(evt,pb)] = -st
                    self.conf_pbs[(evt,pb)] = -en + 1

            if clust > self.max_clust:
                self.max_clust = clust

        self.mapped = conf_clust is not None

        #Handle unmapped reads
        #TODO: mark them more explicitly
        if not self.mapped:
                #conf_clust = self.max_clust
                #self.conf_evt = len(self.events)-1
            #else:
            if self.max_clust is None:
                self.conf_evt = self.max_evt

        self._set_conf_clust(conf_clust)
        self._parse_expired(clust_ids)

        self.seeds_loaded = True
        return True

    #sets confident cluster related vars
    def _set_conf_clust(self, cc):
        self.conf_clust = cc

        if cc is None: 
            self.plot_conf = False
            return


        self.conf_idx = np.searchsorted(cc.evts, self.conf_evt, side='right')-1
        self.conf_len = cc.lens[self.conf_idx]

        #TODO use cc.evrf_ratio(e,r)
        evt_span = self.conf_evt - cc.evt_st

        if self.max_evt is None:
            max_max_evt = self.conf_evt + (evt_span * self.CONF_PAD_COEF)
            self.max_evt = int(np.round(min(cc.evt_en+1, max_max_evt)))

        #TODO rename min/max_clust_ref
        self.min_idx = np.searchsorted(cc.evts, self.min_evt, side='left')
        self.max_idx = np.searchsorted(cc.evts, self.max_evt, side='right')-1

        rst = cc.rsts[self.min_idx] - SEED_LEN
        ren = cc.rens[self.max_idx]+1

        #ref_span = ren-rst
        #if ref_span > 0:
        #    self.evrf_ratio = (self.max_evt) / ref_span
        #else:
        #    self.evrf_ratio = 1

        self.max_clust_len = cc.lens[self.max_idx]

        if self.mm2_paf is not None:
            mst,men = self.mm2_paf.ext_ref()
            self.plot_conf = self.ref_name == cc.rf
            if coords_inverted(cc.fwd):
                self.plot_conf = self.plot_conf and max(ren,mst) < min(rst,men)
            else:
                self.plot_conf = self.plot_conf and max(rst,mst) < min(ren,men)
        else:
            self.plot_conf = True
            self.fwd = cc.fwd

            if coords_inverted(cc.fwd):
                self.min_ref = -rst
                self.max_ref = -ren
            else:
                self.min_ref = rst
                self.max_ref = ren

            self.ref_name = cc.rf
            self.inv_coords = coords_inverted(cc.fwd)
            #if self.min_ref is None:
            #if self.max_ref is None:
            #if self.ref_name is None:


    #finds clusters that should expire
    #by the end of loaded events
    def _parse_expired(self, clust_ids):
        for cid in clust_ids:
            c = self.clusts[cid]
            if c.expired(self.max_evt):
                self.clusts_exp.append(c)
                del self.clusts[cid]

    def _parse_moves(self, moves):
        return np.flip()


    def parse_paths(self):
        if self.paths_loaded: return False

        path_counts = list()

        paths = open(self.paths_fname)
        head_tabs = paths.readline().split()
        C = {head_tabs[i] : i for i in range(len(head_tabs))}

        for line in paths:
            tabs = line.split()
            if tabs[0] == head_tabs[0]: continue
            path_id  = tuple(map(int, tabs[C['id']].split(':')))

            ref_st = self.conf_pbs.get(path_id, None)

            is_conf = ref_st is not None

            if not is_conf and self.max_path_fm == 0: continue

            evt, pb = path_id
            moves = [bool(c=='1') for c in tabs[C['moves']]]

            if evt < self.min_evt or evt >= self.max_evt: 
                continue

            evt = evt-self.min_evt

            if not self.evt_mask[evt-len(moves)]: 
                continue

            evts = np.arange(evt-len(moves), evt) + 1

            if not is_conf:
                fm_len    =   int(tabs[C['fm_len']])

                if fm_len > self.max_path_fm: continue

                fm_start  =   int(tabs[C['fm_start']])
                full_len  =   int(tabs[C['full_len']])

                i = np.searchsorted(self.range_fms, fm_start)

                while i < len(self.range_fms) and self.range_fms[i] < fm_start+fm_len: 
                    #if True:#not self.bc_loaded or evt in self.range_fm_evts[i]:
                    if not self.bc_loaded or evt in self.range_fm_evts[i]:
                        ref_en = self.fm_to_ref[self.range_fms[i]] - 3

                        #print("%d:%d\t%d\t%d" % (evt,pb,ref_en,fm_len))
                        print("%d:%d" % (evt,pb))

                        #refs = ref_en - np.sum(moves) + np.cumsum(np.flip(moves)) - 3

                        p = (evt,ref_en)
                        self.path_lens[(evt,ref_en)] = min(
                                np.log2(fm_len), self.path_lens.get(p, np.inf)
                        )

                        #self.path_lens[(evt,ref_en)] = min(self.path_lens[(evt,ref_en)], np.log2(fm_len))

                        #for e,r in zip(evts,refs):
                        #    if e < 0: continue
                        #    self.path_lens[(e,r)] = min(self.path_lens[(e,r)], np.log2(fm_len))
                    i += 1

            else:
                refs = ref_st + (np.cumsum(np.flip(moves)) - 1)

                for e,r in zip(evts,refs):
                    if e < 0: continue
                    self.conf_dots.add( (e,r) )

            #kmer       =       tabs[C['kmer']]
            #match_prob = float(tabs[C['match_prob']])
            #self.seed_kmers[(e, refs[-1])] = (kmer, match_prob)


            #print(path_id,parent,fm_start,fm_len,kmer,full_len,seed_prob,moves)

        #    #Store number of paths at each event position
        #    while len(path_counts) - 1 < evt:
        #        path_counts.append(0)
        #    path_counts[-1] += 1

        #self.paths_loaded = True

        return True

    def parse_conf(self):
        if self.conf_loaded: return False

        self.conf_evts = list()
        self.conf_clusts = list()
        self.conf_tops = list()
        self.conf_means = list()

        conf_in = open(self.conf_fname)
        conf_in.readline()

        for line in conf_in:
            evt,clust,top,mean = line.split()
            self.conf_evts.append(int(evt))
            self.conf_clusts.append(int(clust))
            self.conf_tops.append(float(top))
            self.conf_means.append(float(mean))
            

    def parse_bc_aln(self, bce_moves):

        #List of ref coords for each query (read) cord
        qr_to_rfs = self._cig_query_to_refs(self.mm2_paf)

        #basecalled event start, end, and stride
        bce_samp_st, bce_moves_pac = bce_moves
            
        bce_samp_en = bce_samp_st + len(bce_moves) * BCE_STRIDE

        bce_st = 0
        bce_en = int((self.max_samp-bce_samp_st+1) // BCE_STRIDE)
        bce_moves = np.unpackbits(bce_moves_pac)[bce_st:bce_en]

        #Read coord of each basecalled event
        bce_qrs = np.cumsum(bce_moves)
        i = np.searchsorted(bce_qrs, self.mm2_paf.qr_st)
        
        #bce_evts = list()
        bce_samps = list()
        bce_refs = list()

        samp = bce_samp_st
        for qr in bce_qrs:
            if samp >= self.min_samp:
                for rf in qr_to_rfs[qr]:
                    bce_samps.append(samp)
                    bce_refs.append(rf)

            samp += BCE_STRIDE

        self.bce_samps = np.array(bce_samps)
        self.bce_refs = np.array(bce_refs) - BCE_K + 1
        self.max_ref = self.min_ref + max(bce_refs)
        
        self.bc_loaded = True


    def _cig_query_to_refs(self, paf):
        cig = paf.tags.get('cg', (None,)*2)[0]
        if cig is None: return None

        qr_rfs = defaultdict(list)

        if IS_RNA:
            qr_i = paf.qr_len - paf.qr_en 
        else:
            qr_i = paf.qr_st

        rf_i = 0#paf.rf_st

        cig_ops = CIG_RE.findall(cig)

        if coords_inverted(paf.is_fwd):
            cig_ops = reversed(cig_ops)

        for l,c in cig_ops:
            l = int(l)
            incr_qr = c in CIG_INCR_RD
            incr_rf = c in CIG_INCR_RF
            qr_j = qr_i + (l if incr_qr else 1)
            rf_j = rf_i + (l if incr_rf else 1)
            #if c == "M":
            for qr,rf in zip(range(qr_i, qr_j), range(rf_i, rf_j)):
                qr_rfs[qr].append(rf)
            #for qr,rf in zip(range(qr_i, qr_j), range(rf_i, rf_j)):
            #    qr_rfs[qr].append(rf)

            if incr_qr:
                qr_i = qr_j 

            if incr_rf:
                rf_i = rf_j 

        return qr_rfs

class SeedCluster:
    EXPIRE_COEF = 5.5

    #TODO add start evt at st-SEED_LEN?
    #would make reference len easier to breakdown
    #and could be nice for plotting?
    def __init__(self, evt, rf, st, en, fwd, cid):
        self.id = cid
        self.rf = rf
        self.evts = [evt]
        self.blocks = [(st, en)] #TODO only use rsts?
        self.gains = [en-st]

        self.rsts = [st] #TODO consider strand
        self.rens = [en] #TODO consider strand

        self.lens = [self.gains[0]]
        self.fwd = fwd
        self.exp_evt = None

        #if IS_RNA:
        #    self.fwd = not fwd


    def expired(self, evt=None):
        return False #

        #if evt == None:
        #    return self.exp_evt is not None

        #if evt - self.evt_en > self.ref_len * self.EXPIRE_COEF:
        #    if self.exp_evt is None:
        #        self.exp_evt = int(np.round(
        #            self.evt_en + (self.ref_len * self.EXPIRE_COEF)
        #        ))
        #    return True
        #return False

    def add_gain(self, gain):
        self.gains.append(gain)
        self.lens.append(len(self)+max(gain,0))

    def add_seed(self, evt, rf, st, en):

        bst, ben = self.blocks[-1]

        self.evts.append(evt)
        self.rsts.append(st)
        self.rens.append(en)

        if (self.rf != rf or max(bst, st) > min(ben, en)):
            self.blocks.append( (st,en) )
            self.add_gain(en-st)

        else:
            l1 = ben-bst
            self.blocks[-1] = (min(st, bst), max(en, ben))
            l2 = self.blocks[-1][1] - self.blocks[-1][0]
            self.add_gain(l2-l1)

        return True

    #TODO add start evt at st-SEED_LEN?
    @property
    def ref_len(self):
        return self.lens[-1]

    @property
    def evt_st(self):
        return self.evts[0]-SEED_LEN

    @property
    def evt_en(self):
        return self.evts[-1]

    def evt_len(self, i=0, j=-1):
        return self.evts[j] - self.evts[i] + 1

    @property
    def evrf_ratio(self, i=0, j=-1):
        return (self.evts[j] - self.evt_st+1) / self.lens[i]

    def __str__(self):
        MAX_REF = 22
        if len(self.rf) < MAX_REF:
            ref = self.rf
        else:
            ref = self.rf[:MAX_REF] + "."

        return "%s:%d-%d (%s)" % (
                ref,
                self.blocks[0][0],
                self.blocks[-1][1],
                "+" if self.fwd else "-")


    #TODO bind confidence to len
    #then make it more complicated
    def __gt__(self, c):
        return c is None or len(self) > len(c)

    def __eq__(self, rhs):
        return self.id == getattr(rhs, "id", None)

    def __len__(self):
        return self.lens[-1]
