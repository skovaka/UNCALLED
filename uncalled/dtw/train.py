from . import Tracks
from .dtw import dtw_pool_iter
import sys
from time import time

def train(conf):
    conf.fast5_reader.load_bc = True
    conf.tracks.load_fast5s = True
    conf.bc_cmp = True
    print("Loading tracks")

    tracks = Tracks(conf=conf)
    
    trainer = tracks.output

    print("setting")
    if conf.train.append:
        conf.pore_model.name = trainer.model.name
    else:
        trainer.set_model(tracks.model)

    if conf.train.skip_dtw:
        print("skipping")
        model_file = trainer.next_model(True)
        return
        #itr = range(conf.train.iterations, conf.train.iterations+1)

    t = time()
    for i in range(conf.train.iterations):
        i = 0
        for chunk in dtw_pool_iter(tracks):
            i += len(chunk)
            trainer.write_buffer(chunk)
            if trainer.is_full():
                break
            print("reads", i, time()-t)
            sys.stdout.flush()
        print("finish", trainer.is_full())

        model_file = trainer.next_model()
        tracks.bam_in.reset()
        conf.pore_model.name = model_file
        tracks = Tracks(conf=conf)

    tracks.close()

