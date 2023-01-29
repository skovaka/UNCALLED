from . import Tracks
from .dtw import dtw_pool_iter

def train(conf):
    conf.fast5_reader.load_bc = True
    conf.tracks.load_fast5s = True
    conf.bc_cmp = True

    tracks = Tracks(conf=conf)
    
    trainer = tracks.output
    trainer.set_model(tracks.model)

    for i in range(conf.train.iterations):
        print("ITER", i)
        for chunk in dtw_pool_iter(tracks):
            trainer.write_buffer(chunk)

            if trainer.is_full():
                break

        model_file = trainer.next_model()
        tracks.bam_in.reset()
        conf.pore_model.name = model_file
        tracks = Tracks(conf=conf)

    tracks.close()

