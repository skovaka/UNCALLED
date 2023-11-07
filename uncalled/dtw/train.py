from . import Tracks
from .dtw import DtwPool, GuidedDTW
import sys
from time import time
from _uncalled import EventDetector
import numpy as np
import pandas as pd
from collections import Counter
import multiprocessing as mp


from ..pore_model import PoreModel, PoreModelParams

def init_model(tracks, k):
    p = tracks.conf.event_detector
    p.min_mean = 0
    p.max_mean = 1000000
    evdt = EventDetector(p)

    currents = list()
    length = 0

    for sam in tracks.bam_in.iter_sam():

        mv = np.array(sam.get_tag("mv"))
        mv_stride = mv[0]
        st = sam.get_tag("ts")
        moves = mv[1:]

        en = st + (np.sum(moves) * mv_stride)

        read = tracks.read_index[sam.query_name]

        c = evdt.get_means(read.signal.to_numpy()[st:en])

        currents.append(c)
        length += len(c)
        if length >= tracks.conf.train.init_events:
            break
    currents = np.concatenate(currents)

    mn = np.mean(currents)
    sd = np.std(currents)
    coef = 3
    cmin = mn-sd*coef
    cmax = mn+sd*coef

    currents = currents[(currents >= cmin) & (currents <= cmax)]

    mn = np.mean(currents)
    sd = np.std(currents)
    coef = 3
    cmin = mn-sd*coef
    cmax = mn+sd*coef

    tracks.conf.pore_model.pa_mean = np.mean(currents)
    tracks.conf.pore_model.pa_stdv = np.std(currents)
    #tracks.set_model(PoreModel((tracks.conf.pore_model, tracks.conf.normalizer.tgt_mean, tracks.conf.normalizer.tgt_stdv)))
    model = PoreModel(params=tracks.conf.pore_model)
    print(model.K, "INIITNTITITITs")
    tracks.set_model(model)

def train(conf):
    mp.set_start_method("spawn")
    #conf.tracks.load_signal = True
    conf.tracks.layers.append("moves")
    conf.mvcmp = True

    prms = conf.train

    if len(prms.init_model) > 0:
        conf.pore_model.name = prms.init_model
        #conf.pore_model.k = prms.
    elif not prms.append and prms.kmer_len is not None:
        conf.pore_model.k = prms.kmer_len
        orig_dtw_iters = conf.dtw.iterations
        conf.dtw.iterations = 0
    elif not prms.append:
        raise ValueError(f"Must define kmer_length, init_model, or run in append mode")

    #tracks = Tracks(model=PoreModel(prms.init_model), conf=conf)
    tracks = Tracks(model=PoreModel(params=conf.pore_model), conf=conf)
    _ = tracks.read_index.default_model

    if conf.dtw.iterations == 0:
        init_model(tracks, prms.kmer_len)

    trainer = tracks.output
    trainer.model = tracks.model

    if prms.append:
        conf.pore_model.name = trainer.model.name
    else:
        trainer.set_model(tracks.model)

    print(conf.train.init_mode)

    if prms.skip_dtw:
        model_file = trainer.next_model(True)
        return
        #itr = range(prms.iterations, prms.iterations+1)

    bam_in = tracks.bam_in.input
    bam_in.reset()
    
    t = time()
    for i in range(prms.iterations):
        print("iter", i+1, "of", prms.iterations, tracks.index.model.K)
        bam_start = bam_in.tell()

        #if i > 2:
        #    tracks.conf.dtw.skip_cost = 2
        #    tracks.conf.tracks.mask_skips = None
        #    tracks.conf.dtw.iterations = 2

        status_counts = Counter()
        #for aln in tracks.bam_in.iter_alns():
        #    dtw = GuidedDTW(tracks, aln)
        #    status_counts[dtw.status] += 1
        pool = DtwPool(tracks)
        for chunk,counts in pool: #dtw_pool_iter(tracks):
            status_counts.update(counts)
            print("LOOP")
            trainer.write_buffer(chunk)

            if trainer.is_full():
                break

        print("hereo")

        pool.close()

        #If EOF reached, reset BAM and read until full or entire file was read
        if not trainer.is_full():
            bam_in.reset()

            #for aln in tracks.bam_in.iter_alns():
            #    dtw = GuidedDTW(tracks, aln)
            #    status_counts[dtw.status] += 1
            pool = DtwPool(tracks)
            for chunk,counts in pool: #dtw_pool_iter(tracks):
                status_counts.update(counts)
                trainer.write_buffer(chunk)
                if bam_in.tell() >= bam_start or trainer.is_full():
                    break
            pool.close()

        prms.append = False
        model = trainer.next_model()
        if tracks.conf.dtw.iterations == 0:
            tracks.conf.dtw.iterations = orig_dtw_iters

        tracks.set_model(model)

        sys.stderr.write(str(status_counts)+"\n")

    tracks.close()

