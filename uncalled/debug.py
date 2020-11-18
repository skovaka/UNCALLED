#!/usr/bin/env python

import sys, os
import numpy as np
import argparse
import bisect
from scipy.stats import linregress
from file_read_backwards import FileReadBackwards
from collections import defaultdict

SEED_LEN = 22
MIN_CLUST = 25
SAMPLE_RATE = 4000
CHUNK_LEN = 4000

MAX_CHUNK_DEF = 3

class DebugParser:
    CONF_PAD_COEF = 2

    def __init__(self, 
                 prefix, 
                 paf, 
                 min_chunk=None,
                 max_chunk=None,
                 min_samp=None,
                 max_samp=None,
                 tgt_cid=None,
                 load_seeds=True,
                 load_events=True,
                 load_paths=True,
                 load_conf=True):

        self.rid = paf.qr_name

        #Cofident seed cluster ID, confident event, and reference length
        #(conf_evt is where mapping would normally end)
        self.conf_cid = paf.tags.get('sc', (None,)*2)[0]
        self.conf_evt = paf.tags.get('ce', (None,)*2)[0]
        self.conf_samp = None

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

        print(self.max_chunk, max_chunk, max_samp)

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
        self.min_clust_ref = None #TODO name span?
        self.max_clust_ref = None #TODO name span?

        #TODO clean these up
        self.conf_pbs = dict()
        self.dots = defaultdict(int)
        self.seed_kmers = dict()
        self.conf_ref_bounds = None

        self.evts_loaded = False
        if load_events:
            self.evts_in = open(prefix + self.rid + "_events.tsv")
            self.parse_events()

        self.seeds_loaded = False
        if load_seeds:
            self.seeds_in  = open(prefix + self.rid + "_seeds.bed")
            self.parse_seeds()

        #TODO: paths
        self.paths_loaded = False
        if load_paths:
            self.paths_fname = prefix + self.rid + "_paths.tsv"
            self.parse_paths()

        self.conf_loaded = False
        if load_conf:
            self.conf_fname = prefix + self.rid + "_conf.tsv"
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
                print(en, self.min_samp)
                self.min_evt = evt
                print(self.min_evt, "MIN")

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

        if self.max_evt is None or self.max_evt >= evt: 
            self.max_evt = evt

        if self.max_samp is None:
            self.max_samp = st+ln
            self.max_chunk = (self.max_samp-1) // CHUNK_LEN

        self.evts_loaded = True

        self.win_means = np.array(self.win_means)
        self.win_stdvs = np.array(self.win_stdvs)
        self.evt_mask = np.array(self.evt_mask)
        self.evt_mask_map = np.array(self.evt_mask_map)

        return True

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
                if clust.fwd:
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
            if self.max_clust is not None:
                conf_clust = self.max_clust
                self.conf_evt = conf_clust.evts[-1]
            else:
                self.conf_evt = self.max_evt

        self._set_conf_clust(conf_clust)
        self._parse_expired(clust_ids)

        self.seeds_loaded = True
        return True

    #sets confident cluster related vars
    def _set_conf_clust(self, cc):
        self.conf_clust = cc

        if cc is None: return

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

        if cc.fwd:
            rst = cc.rsts[self.min_idx] - SEED_LEN
            ren = cc.rens[self.max_idx]+1
        else:
            rst = -cc.rsts[self.min_idx] - SEED_LEN
            ren = -cc.rens[self.max_idx]-1

        ref_span = ren-rst
        if ref_span > 0:
            self.evrf_ratio = (self.max_evt) / ref_span
        else:
            self.evrf_ratio = 1

        self.max_clust_len = cc.lens[self.max_idx]

        self.min_clust_ref = rst
        self.max_clust_ref = rst + (self.max_evt / self.evrf_ratio)

    #finds clusters that should expire
    #by the end of loaded events
    def _parse_expired(self, clust_ids):
        for cid in clust_ids:
            c = self.clusts[cid]
            if c.expired(self.max_evt):
                self.clusts_exp.append(c)
                del self.clusts[cid]

    def parse_paths(self):
        if self.paths_loaded: return False

        path_counts = list()

        paths_fwd = open(self.paths_fname)
        head_tabs = paths_fwd.readline().split()
        C = {head_tabs[i] : i for i in range(len(head_tabs))}

        #paths_fwd.close()
        #paths_rev = FileReadBackwards(self.paths_fname)
        #for line in paths_rev:

        for line in paths_fwd:
            tabs = line.split()
            if tabs[0] == head_tabs[0]: continue
            path_id  = tuple(map(int, tabs[C['id']].split(':')))

            ref_st = self.conf_pbs.get(path_id, None)
            if ref_st is None: continue

            evt, pb = path_id

            moves = [bool(c=='1') for c in tabs[C['moves']]]

            if evt < self.min_evt or evt >= self.max_evt: 
                continue

            e = evt-self.min_evt

            if not self.evt_mask[e-len(moves)]: 
                continue

            evts = np.arange(e-len(moves), e) + 1
            refs = ref_st + np.cumsum(np.flip(moves))

            for e,r in zip(evts,refs):
                if e >= 0: 
                    self.dots[(e,r)] = True

            kmer      =       tabs[C['kmer']]
            self.seed_kmers[(evt, refs[-1])] = kmer

            #fm_start  =   int(tabs[C['fm_start']])
            #fm_len    =   int(tabs[C['fm_len']])
            #full_len  =   int(tabs[C['full_len']])
            #seed_prob = float(tabs[C['seed_prob']])

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
