#!/usr/bin/env python

import sys, os
import numpy as np
import argparse
import bisect
from scipy.stats import linregress

SEED_LEN = 22
MIN_CLUST = 25
SAMPLE_RATE = 4000
CHUNK_LEN = 4000

class DebugParser:
    CONF_PAD_COEF = 2

    def __init__(self, 
                 prefix, 
                 paf, 
                 chunk_start=0,
                 chunk_count=None,
                 load_seeds=True,
                 load_events=True,
                 load_paths=True):

        self.rid = paf.qr_name

        #Cofident seed cluster ID, confident event, and reference length
        #(conf_evt is where mapping would normally end)
        self.conf_cid = paf.tags.get('sc', (None,)*2)[0]
        self.conf_evt = paf.tags.get('ce', (None,)*2)[0]

        self.chunk_st = chunk_start
        self.samp_st = self.chunk_st * SAMPLE_RATE

        #Set bounds based on argument
        if chunk_count is not None:
            self.chunk_en = chunk_start + chunk_count
            self.samp_en = self.chunk_en * SAMPLE_RATE

        #Set to constant if unmapped
        elif self.conf_evt is None:
            self.chunk_en = self.chunk_st + 3 #TODO store const
            self.samp_en = self.chunk_en * SAMPLE_RATE

        #If mapped, TBD in parse evts
        else:
            self.chunk_en = self.samp_en = None

        self.chunk_evt_bounds = [0]

        self.min_evt = 0 if chunk_start == 0 else None
        self.max_evt = None

        #TBD in parse_seeds
        self.conf_len = None 
        self.conf_clust = None
        self.max_clust = None

        self.max_clust_len = None
        self.max_clust_ref = None #TODO name span?

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
            self.paths_in  = open(prefix + self.rid + "_paths.tsv")
            self.parse_paths()

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

            if st >= next_chunk_samp:
                chunk += 1
                next_chunk_samp += CHUNK_LEN
                self.chunk_evt_bounds.append(evt)

            if chunk < self.chunk_st: continue
            if self.chunk_en is not None and chunk >= self.chunk_en: break
            
            self.events.append( (st,ln,mn,sd) )
            self.norms.append( (norm_sc,norm_sh) )

            self.win_means.append(win_mn)
            self.win_stdvs.append(win_sd)
            self.evt_mask.append(mask == 1)

            if self.evt_mask[-1]:
                self.evt_mask_map.append(evt)

            if self.chunk_en is None and self.conf_evt == evt:
                self.chunk_en = chunk+1

            evt += 1

            if self.samp_en:
                break

        if self.max_evt is None or self.max_evt >= evt: 
            self.max_evt = evt
            self.samp_en = st+ln

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

            if clust > self.max_clust:
                self.max_clust = clust

        #Handle unmapped reads
        #TODO: mark them more explicitly
        if conf_clust == None:
            conf_clust = self.max_clust
            self.conf_evt = conf_clust.evts[-1]

        self._set_conf_clust(conf_clust)
        self._parse_expired(clust_ids)

        self.seeds_loaded = True
        return True

    #sets confident cluster related vars
    def _set_conf_clust(self, cc):
        self.conf_clust = cc
        self.conf_idx = np.searchsorted(cc.evts, self.conf_evt, side='right')-1
        self.conf_len = cc.lens[self.conf_idx]

        #TODO use cc.evrf_ratio(e,r)
        evt_span = self.conf_evt - cc.evt_st

        if self.max_evt is None:
            max_max_evt = self.conf_evt + (evt_span * self.CONF_PAD_COEF)
            self.max_evt = int(np.round(min(cc.evt_en+1, max_max_evt)))

        self.max_idx = np.searchsorted(cc.evts, self.max_evt, side='right')-1

        if cc.fwd:
            rst = cc.rsts[0]
            ren = cc.rens[self.max_idx]+1
        else:
            rst = cc.rsts[self.max_idx]
            ren = cc.rens[0]+1

        ref_span = ren-rst
        self.evrf_ratio = evt_span / ref_span

        #self.max_clust_len = int(np.round(self.max_evt / self.evrf_ratio))
        self.max_clust_len = cc.lens[self.max_idx]

        self.min_ref = rst
        self.max_clust_ref = rst+int(np.round(self.max_evt / self.evrf_ratio))

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

        head_tabs = self.paths_in.readline().split()
        C = {head_tabs[i] : i for i in range(len(head_tabs))}

        for line in self.paths_in:
            tabs = line.split()
            path_id   =       tabs[C['id']]
            parent    =       tabs[C['parent']]
            fm_start  =   int(tabs[C['fm_start']])
            fm_len    =   int(tabs[C['fm_len']])
            kmer      =       tabs[C['kmer']]
            full_len  =   int(tabs[C['full_len']])
            seed_prob = float(tabs[C['seed_prob']])
            moves     =       tabs[C['moves']]
            evt, buf_id = map(int, path_id.split(':'))

            #print(path_id,parent,fm_start,fm_len,kmer,full_len,seed_prob,moves)

            #Store number of paths at each event position
            while len(path_counts) - 1 < evt:
                path_counts.append(0)
            path_counts[-1] += 1

        self.paths_loaded = True

        return True

class SeedCluster:
    EXPIRE_COEF = 5.5

    #TODO add start evt at st-SEED_LEN?
    #would make reference len easier to breakdown
    #and could be nice for plotting?
    def __init__(self, evt, rf, st, en, fwd, cid):
        self.id = cid
        self.rf = rf
        self.evts = [evt]
        self.blocks = [(st, en)]
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
        return "%s\t%d\t%d\t%d\t%d" % (
                self.rf,
                self.blocks[0][0],
                self.blocks[-1][1],
                self.id,
                len(self))


    #TODO bind confidence to len
    #then make it more complicated
    def __gt__(self, c):
        return c is None or len(self) > len(c)

    def __eq__(self, rhs):
        return self.id == getattr(rhs, "id", None)

    def __len__(self):
        return self.lens[-1]
