from . import Tracks
from .dtw import dtw_pool_iter
import sys
from time import time

from ..pore_model import PoreModel, PoreModelParams

def train(conf):
    conf.fast5_reader.load_bc = True
    conf.tracks.load_fast5s = True
    conf.bc_cmp = True

    prms = conf.train

    if prms.init_model is not None:
        conf.pore_model.name = prms.init_model
    elif prms.kmer_len is not None:
        #conf.pore_model.name = ""
        conf.pore_model.k = prms.kmer_len
        conf.dtw.iterations = 0
    elif not prms.append:
        raise ValueError(f"Must define kmer_length, init_model, or run in append mode")

    tracks = Tracks(conf=conf)
    
    trainer = tracks.output

    if prms.append:
        conf.pore_model.name = trainer.model.name
    else:
        trainer.set_model(tracks.model)

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
            print("BAM RESET", bam_start)
            bam_in.reset()
            for chunk in dtw_pool_iter(tracks):
                count += len(chunk)
                trainer.write_buffer(chunk)
                if bam_in.tell() >= bam_start or trainer.is_full():
                    print("done", trainer.is_full())
                    break

        prms.append = False
        tracks.conf.dtw.iterations = 1

        tracks.set_model(trainer.next_model())

    tracks.close()

