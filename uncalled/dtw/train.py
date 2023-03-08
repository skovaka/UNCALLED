from . import Tracks
from .dtw import dtw_pool_iter
import sys
from time import time
from _uncalled import EventDetector
import numpy as np

from ..pore_model import PoreModel, PoreModelParams

def init_model(tracks, k):
    p = tracks.conf.event_detector
    p.min_mean = 0
    p.max_mean = 250
    evdt = EventDetector(p)

    currents = list()
    length = 0

    for read in tracks.read_index:
        st = read.template_start
        en = st + np.sum(read.moves)*5
        c = evdt.get_means(read.signal[st:en])
        currents.append(read.signal[st:en])
        length += len(c)
        if length >= tracks.conf.train.init_events:
            break
    currents = np.concatenate(currents)

    tracks.conf.normalizer.tgt_mean = np.mean(currents)
    tracks.conf.normalizer.tgt_stdv = np.std(currents)
    tracks.set_model(PoreModel((tracks.conf.pore_model, tracks.conf.normalizer.tgt_mean, tracks.conf.normalizer.tgt_stdv)))

def train(conf):
    conf.tracks.load_fast5s = True
    conf.mvcmp = True

    prms = conf.train

    if prms.init_model is not None:
        conf.pore_model.name = prms.init_model
    elif prms.kmer_len is not None:
        conf.pore_model.name = ""
        conf.pore_model.k = prms.kmer_len
        conf.pore_model.shift = (prms.kmer_len-1) // 2
        print(conf.pore_model.k)
        conf.dtw.iterations = 0
    elif not prms.append:
        raise ValueError(f"Must define kmer_length, init_model, or run in append mode")

    tracks = Tracks(conf=conf)

    if conf.dtw.iterations == 0:
        init_model(tracks, prms.kmer_len)

    trainer = tracks.output
    trainer.model = tracks.model

    if prms.append:
        conf.pore_model.name = trainer.model.name
    else:
        trainer.set_model(tracks.model)

    print(tracks.conf.pore_model.name)

    if prms.skip_dtw:
        model_file = trainer.next_model(True)
        return
        #itr = range(prms.iterations, prms.iterations+1)

    bam_in = tracks.bam_in.input
    
    t = time()
    for i in range(prms.iterations):
        print("iter", i+1, "of", prms.iterations)
        count = 0
        bam_start = bam_in.tell()
        for chunk in dtw_pool_iter(tracks):
            count += len(chunk)
            trainer.write_buffer(chunk)
            if trainer.is_full():
                break

        #If EOF reached, reset BAM and read until full or entire file was read
        if not trainer.is_full():
            bam_in.reset()
            for chunk in dtw_pool_iter(tracks):
                count += len(chunk)
                trainer.write_buffer(chunk)
                if bam_in.tell() >= bam_start or trainer.is_full():
                    break

        prms.append = False
        tracks.conf.dtw.iterations = 1

        tracks.set_model(trainer.next_model())

    tracks.close()

