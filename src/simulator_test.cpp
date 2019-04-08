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

    //std::ifstream reads_file(reads_fname);
    //if (!reads_file) {
    //    std::cerr << "Error: couldn't open '" << reads_fname << "'\n";
    //    return 1;
    //}
    //std::string read_fname;
    //std::vector<std::string> fast5_names;
    //fast5_names.reserve(read_count);
    //u32 i = 0;
    //while (getline(reads_file, read_fname) && i < read_count) {
    //    if (read_fname[0] == '#') continue;
    //    fast5_names.push_back(read_fname);
    //    i++;
    //}
    //std::cerr << t.get() << "\tfilenames read, splitting into chunks...\n";
    
    std::cerr << "Loading reads\n";
    std::cerr.flush();

    std::vector<Fast5Read> reads;
    load_multi_fast5(reads_fname, reads);

    std::cerr << "Reads loaded (" << (t.lap() / 1000) << " sec), sorting\n";
    std::cerr.flush();

    std::sort(reads.begin(), reads.end());
    while (reads.size() > read_count) reads.pop_back();

    std::cerr << "Reads sorted (" << (t.lap() / 1000) << " sec), splitting into chunks\n";
    std::cerr.flush();

    ChunkSim sim(read_count, 512, 4000, speed);
    sim.add_reads(reads);

    std::cerr << "Chunks split (" << (t.lap() / 1000) << " sec), aligning\n";
    std::cerr.flush();

    ChunkPool pool(params, 512, nthreads);

    const u64 MAX_SLEEP = 100;

    while (sim.is_running || !pool.all_finished()) {
        u64 t0 = t.get();
        
        //if (!chunks.empty()) {
        //    std::cout << t.get() << "\t" << chunks.size() << " chunks\n";
        //}

        for (ReadLoc aln : pool.update()) {
            std::cout << aln.str() << "\n";
            std::cout.flush();
            sim.unblock(aln.get_channel(), aln.get_number());
            //sim.stop_receiving_read(aln.get_channel(), aln.get_number());
        }

        std::vector<ChChunk> chunks = sim.get_read_chunks();
        for (ChChunk &ch : chunks) {
            //std::cout << "# chunk " << ch.first << " " << ch.second.id << "\n";
            std::cout.flush();

            u16 channel = ch.first, number = ch.second.number;

            bool last = ch.second.raw_data.size() < 4000;

            pool.add_chunk(ch.first, ch.second);
            //if (!pool.add_chunk(ch.first, ch.second)) {
                //`std::cout << "# couldn't add\n";
            //}
            
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
