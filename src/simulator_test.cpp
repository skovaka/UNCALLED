#include <iostream>
#include <unistd.h>
#include "fast5_reader.hpp"
#include "chunk_pool.hpp"
#include "mapper.hpp"

const std::string MODEL = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers.txt";
const std::string PROBFN = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";


int main(int argc, char** argv) {
    std::string index(argv[1]), reads_fname(argv[2]);
    u32 read_count = atoi(argv[3]);

    MapperParams params(index, MODEL,
                        22,    //seed_len
                        25,    //min_aln_len
                        0,     //min_rep_len
                        50,    //max_rep_copy
                        8,     //max_consec_stay
                        10000, //max_paths
                        30000, //max_events_proc
                        6000,     //evt_buffer_len TODO: what should it be
                        3,     //evt_winlen1
                        6,     //evt_winlen2
                        5,     //evt_batch_size
                        1.4,   //evt_thresh1
                        9.0,   //evt_thresh2
                        0.2,   //evt_peak_height
                        30,    //evt_min_mean
                        150,   //evt_max_mean
                        0.5,   //max_stay_frac
                        -3.75, //min_seed_prob
                        7.00,  //min_mean_conf
                        2.25); //min_top_conf

    std::ifstream reads_file(reads_fname);

    if (!reads_file) {
        std::cerr << "Error: couldn't open '" << reads_fname << "'\n";
        return 1;
    }

    Timer t;

    std::string read_fname;
    std::vector<std::string> fast5_names;
    fast5_names.reserve(read_count);

    u32 i = 0;
    while (getline(reads_file, read_fname) && i < read_count) {
        if (read_fname[0] == '#') continue;
        fast5_names.push_back(read_fname);
        i++;
    }
    std::cerr << t.get() << "\tfilenames read, splitting into chunks...\n";
    
    ChunkSim sim(read_count, 512, 4000, fast5_names);
    ChunkPool pool(params, 512, 1);

    std::cerr << t.get() << "\tchunks generated, aligning\n";

    const u64 MAX_SLEEP = 250;

    while (sim.is_running) {
        u64 t0 = t.get();
        std::vector<ChChunk> chunks = sim.get_read_chunks();
        
        //if (!chunks.empty()) {
        //    std::cout << t.get() << "\t" << chunks.size() << " chunks\n";
        //}

        for (ChChunk &ch : chunks) {
            std::cout << "# chunk " << ch.first << " " << ch.second.id << "\n";
            std::cout.flush();
            if (!pool.add_chunk(ch.first, ch.second)) {
                std::cout << "Couldn't add\n";
            }
        }


        for (ReadLoc aln : pool.update()) {
            std::cout << aln.str() << "\n";
            std::cout.flush();
            sim.unblock(aln.get_channel(), aln.get_number());
        }

        u64 dt = t.get() - t0;
        if (dt < MAX_SLEEP) usleep(1000*(MAX_SLEEP - dt));
    }

}
