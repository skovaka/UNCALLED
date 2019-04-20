#include <iostream>
#include <unistd.h>
#include "fast5_reader.hpp"
#include "chunk_pool.hpp"
#include "mapper.hpp"

//const std::string MODEL = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers.txt";
//const std::string PROBFN = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

const std::string MODEL = "/home-4/skovaka1@jhu.edu/code/UNCALLED/src/uncalled/models/r94_5mers.txt";
const std::string PROBFN = "/home-4/skovaka1@jhu.edu/code/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

int main(int argc, char** argv) {
    std::string index(argv[1]), reads_fname(argv[2]);
    u32 nthreads = atoi(argv[3]), read_count = atoi(argv[4]);
    float speed = atof(argv[5]), evt_timeout = atof(argv[6]);
    u32 max_chunks_proc = atoi(argv[7]);
    
    MapperParams params(index, MODEL,
                        22,    //seed_len
                        25,    //min_aln_len
                        0,     //min_rep_len
                        50,    //max_rep_copy
                        8,     //max_consec_stay
                        10000, //max_paths
                        30000,
                        max_chunks_proc, //max_events_proc
                        6000,     //evt_buffer_len TODO: what should it be
                        3,     //evt_winlen1
                        6,     //evt_winlen2
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
                        2.25); //min_top_conf

    Timer t;

    //std::cerr << t.get() << "\tfilenames read, splitting into chunks...\n";
    
    //std::cerr << "Loading reads\n";
    //std::cerr.flush();

    //load_multi_fast5(reads_fname, reads);

    //std::cerr << "Reads loaded (" << (t.lap() / 1000) << " sec), sorting\n";
    //std::cerr.flush();

    //std::sort(reads.begin(), reads.end());
    //while (reads.size() > read_count) reads.pop_back();

    //std::cerr << "Reads sorted (" << (t.lap() / 1000) << " sec), splitting into chunks\n";
    //std::cerr.flush();

    std::ifstream reads_file(reads_fname);
    if (!reads_file) {
        std::cerr << "Error: couldn't open '" << reads_fname << "'\n";
        return 1;
    }

    std::string read_fname;
    std::vector<std::string> fast5_names;
    fast5_names.reserve(read_count);
    while (getline(reads_file, read_fname)) {
        if (read_fname[0] == '#') continue;
        fast5_names.push_back(read_fname);
    }

    ChunkSim sim(read_count, 512, 4000, speed);
    sim.add_files(fast5_names);

    std::cerr << "Chunks split (" << (t.lap() / 1000) << " sec), aligning\n";
    std::cerr.flush();

    ChunkPool pool(params, 512, nthreads);

    const u64 MAX_SLEEP = 100;

    sim.start();
    while (sim.is_running() || !pool.all_finished()) {
        u64 t0 = t.get();
        
        //if (!chunks.empty()) {
        //    std::cout << t.get() << "\t" << chunks.size() << " chunks\n";
        //}

        for (ReadLoc aln : pool.update()) {
            if (aln.is_valid()) {
                sim.unblock(aln.get_channel(), aln.get_number());
                aln.set_unblocked();
            } else {
                sim.stop_receiving_read(aln.get_channel(), aln.get_number());
            }
            std::cout << aln.str() << "\n";
            std::cout.flush();
        }

        std::vector<Chunk> chunks = sim.get_read_chunks();
        for (Chunk &ch : chunks) {
            //std::cout << "# chunk " << ch.first << " " << ch.second.id << "\n";
            std::cout.flush();

            u16 channel = ch.get_channel(), number = ch.get_number();

            bool last = ch.size() < 4000;

            //pool.add_chunk(ch);
            if (!pool.add_chunk(ch)) {
                std::cout << "# couldn't add\n";
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
