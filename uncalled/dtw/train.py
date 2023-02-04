from . import Tracks
from .dtw import dtw_pool_iter
import sys
from time import time

def train(conf):
    conf.fast5_reader.load_bc = True
    conf.tracks.load_fast5s = True
    conf.bc_cmp = True

    tracks = Tracks(conf=conf)
    
    trainer = tracks.output

    if conf.train.append:
        conf.pore_model.name = trainer.model.name
    else:
        trainer.set_model(tracks.model)

    if conf.train.skip_dtw:
        model_file = trainer.next_model(True)
        return
        #itr = range(conf.train.iterations, conf.train.iterations+1)

    t = time()
    for i in range(conf.train.iterations):
        print("iter", i+1, "of", conf.train.iterations)
        count = 0
        for chunk in dtw_pool_iter(tracks):
            count += len(chunk)
            trainer.write_buffer(chunk)
            if trainer.is_full():
                print("DONE")
                break
            print("reads", count, time()-t)
            sys.stdout.flush()
        print("finish", trainer.is_full())

        conf.train.append = False

        model_file = trainer.next_model()
        tracks.bam_in.reset()
        conf.pore_model.name = model_file
        tracks = Tracks(conf=conf)

    tracks.close()

