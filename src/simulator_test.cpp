#include <iostream>
#include <unistd.h>
#include "params.hpp"
#include "simulator.hpp"
#include "chunk_pool.hpp"
#include "read_buffer.hpp"

//const std::string MODEL = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers.txt";
//const std::string PROBFN = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

const std::string MODEL = "/home-4/skovaka1@jhu.edu/code/UNCALLED/src/uncalled/models/r94_5mers.txt";
const std::string PROBFN = "/home-4/skovaka1@jhu.edu/code/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

int main(int argc, char** argv) {

    std::string index(argv[1]), reads_fname(argv[2]);
    u32 nthreads = atoi(argv[3]), read_count = atoi(argv[4]);
    float speed = atof(argv[5]), evt_timeout = atof(argv[6]);
    u32 max_chunks_proc = atoi(argv[7]);
    
    Params::init_sim(index, MODEL,
                        22,    //seed_len
                        25,    //min_aln_len
                        0,     //min_rep_len
                        50,    //max_rep_copy
                        8,     //max_consec_stay
                        10000, //max_paths
                        30000, //max_events_proc
                        max_chunks_proc, //max_chunks_proc
                        6000,     //evt_buffer_len TODO: what should it be
                        3,     //evt_winlen1
                        6,     //evt_winlen2
                        nthreads,
                        512,   //channels
                        4000,  //chunklen
                        5,     //evt_batch_size
                        evt_timeout,
                        1.4,   //evt_thresh1
                        9.0,   //evt_thresh2
                        0.2,   //evt_peak_height
                        30,    //evt_min_mean
                        150,   //evt_max_mean
                        0.5,   //max_stay_frac
                        -3.75, //min_seed_prob
                        7.00,  //min_mean_conf
                        2.25,  //min_top_conf
                        speed);

    Timer t;

    Simulator sim;
    sim.add_fast5s(reads_fname, read_count);

    std::cerr << "Chunks split (" << (t.lap() / 1000) << " sec), aligning\n";
    std::cerr.flush();

    ChunkPool pool;

    const u64 MAX_SLEEP = 100;

    sim.start();
    while (sim.is_running() || !pool.all_finished()) {
        u64 t0 = t.get();
        
        //if (!chunks.empty()) {
        //    std::cout << t.get() << "\t" << chunks.size() << " chunks\n";
        //}

        u16 channel;
        u32 number;
        Paf paf;
        for (MapResult m : pool.update()) {
            std::tie(channel, number, paf) = m;
            if (paf.is_mapped()) {
                paf.set_float(Paf::Tag::UNBLOCK, sim.get_time(channel));
                sim.unblock(channel, number);
                //aln.set_unblocked();TODO replace this
            } else {
                sim.stop_receiving_read(channel, number);
                paf.set_float(Paf::Tag::KEEP, sim.get_time(channel));
            }
            paf.print();
        }

        std::vector<Chunk> chunks = sim.get_read_chunks();
        for (Chunk &ch : chunks) {
            //std::cout << "# chunk " << ch.first << " " << ch.second.id << "\n";
            std::cout.flush();

            u16 channel = ch.get_channel(), number = ch.get_number();

            bool last = ch.size() < 4000;

            //pool.add_chunk(ch);
            if (!pool.add_chunk(ch)) {
                std::cerr << "Error: failed to add chunk from " 
                          << ch.get_id() << std::endl;
            }
            
            if (last) { 
                //std::cout << "# ending " << channel << "\n";
                pool.end_read(channel, number);
            }

        }

        u64 dt = t.get() - t0;
        if (dt < MAX_SLEEP) usleep(1000*(MAX_SLEEP - dt));
    }

    pool.stop_all();
    std::cerr << "Reads aligned (" << (t.lap() / 1000) << " sec)\n";
}
