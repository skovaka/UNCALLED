from . import Tracks
from .dtw import dtw_pool_iter

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

    for i in range(conf.train.iterations):
        i = 0
        for chunk in dtw_pool_iter(tracks):
            i += len(chunk)
            print(i, "first")
            trainer.write_buffer(chunk)
            print(i, "next")
            if trainer.is_full():
                print("DONE", i)
                break
        print("finish", trainer.is_full())

        model_file = trainer.next_model()
        tracks.bam_in.reset()
        conf.pore_model.name = model_file
        tracks = Tracks(conf=conf)

    tracks.close()

