#!/usr/bin/env python

import sys, os
import numpy as np
import argparse
import bisect
from scipy.stats import linregress

SEED_LEN = 22
MIN_CLUST = 25
SAMPLE_RATE = 4000

class DebugParser:
    CONF_PAD_COEF = 2

    def __init__(self, prefix, paf, tmplst_sec=None, min_event=0, max_event=None, load_seeds=True, load_events=True, load_paths=True):
        self.rid = paf.qr_name

        #Cofident seed cluster ID, confident event, and reference length
        #(conf_evt is where mapping would normally end)
        self.conf_cid = paf.tags.get('sc', (None,)*2)[0]
        self.conf_evt = paf.tags.get('ce', (None,)*2)[0]

        #TBD in parse_seeds
        self.conf_len = None 
        self.conf_clust = None
        self.max_clust = None

        #Stores maximum bounds to load
        #max_evt/len is set when seeds parsed, if not specified in arg
        self.min_evt = min_event
        self.max_evt = max_event
        self.max_len = None

        #Template start sample
        #maybe should store elsewhere?
        self.tmplst_samp = int(np.round(tmplst_sec*SAMPLE_RATE)) if tmplst_sec is not None else None

        self.seeds_loaded = False
        if load_seeds:
            self.seeds_in  = open(prefix + self.rid + "_seeds.bed")
            self.parse_seeds()

        self.evts_loaded = False
        if load_events:
            self.evts_in = open(prefix + self.rid + "_events.tsv")
            self.parse_events()

        #TODO: paths
        self.paths_loaded = False
        if load_paths:
            self.paths_in  = open(prefix + self.rid + "_paths.tsv")
            self.parse_paths()

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

            clust = self.clusts.get(cid, None)
            replace = (
                clust is not None and 
                clust.expired(evt) and
                (cid != self.conf_cid or evt < self.conf_evt)
            )

            if clust == None or replace:
                if replace:
                    if cid == self.conf_cid:
                        print("Exp: ", evt-clust.evt_en, clust.ref_len)
                    self.clusts_exp.append(clust)

                clust = SeedCluster(evt,rf,st,en,cid)
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
            self.conf_evt = self.conf_clust.evts[-1]

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
        self.evrf_ratio = evt_span / self.conf_len 

        if self.max_evt is None:
            max_max_evt = self.conf_evt + (evt_span * self.CONF_PAD_COEF)
            self.max_evt = int(np.round(min(cc.evt_en, max_max_evt)))

        self.max_len = int(np.round(self.max_evt / self.evrf_ratio))
        self.max_idx = np.searchsorted(cc.evts, self.max_evt, side='right')-1

    #finds clusters that should expire
    #by the end of loaded events
    def _parse_expired(self, clust_ids):
        for cid in clust_ids:
            c = self.clusts[cid]
            if c.expired(self.max_evt):
                self.clusts_exp.append(c)
                del self.clusts[cid]

    def parse_events(self, incl_norm=True):
        if self.evts_loaded: return False
        self.evts_in.readline()

        #Prepare to find template start event
        #TODO too specific to be here?
        self.tmplst_evt = None

        self.events = list()
        self.norms = list()

        self.max_samp = None

        e = 0
        for line in self.evts_in:

            #Skip events if required (min_evt defaults to 0)
            if e < self.min_evt:
                e += 1
                continue
                
            tabs = line.split()
            st,ln = map(int, tabs[:2])
            mn,sd,sc,sh = map(float, tabs[2:])
            
            #Set template start if samp provied
            if (self.tmplst_samp is not None and
                self.tmplst_evt is None and 
                st >= self.tmplst_samp):

                self.tmplst_evt = e

            #Find maximum sample to load
            if (self.max_evt is not None and 
                self.max_samp is None and 
                e >= self.max_evt):

                self.max_samp = st+ln

            self.events.append( (st,ln,mn,sd) )
            self.norms.append( (sc,sh) )

            e += 1

            if (self.tmplst_evt is not None and 
                self.max_samp is not None):

                break

        if self.max_evt is None or self.max_evt >= e: 
            self.max_evt = e
            self.max_samp = st+ln

        self.evts_loaded = True

        return True

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
    def __init__(self, evt, rf, st, en, cid):
        self.id = cid
        self.rf = rf
        self.evts = [evt]
        self.blocks = [(st, en)]
        self.gains = [en-st]
        self.lens = [self.gains[0]]
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
