#include <iostream>
#include <unistd.h>
#include "params.hpp"
#include "simulator.hpp"
#include "chunk_pool.hpp"
#include "read_buffer.hpp"
#include "fast5_pool.hpp"

//const std::string MODEL =  "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers.txt";
//const std::string PROBFN = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

const std::string MODEL = "/home/skovaka/code/UNCALLED/src/uncalled/models/r94_5mers.txt";
const std::string PROBFN = "/home/skovaka/code/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

//const std::string MODEL = "/home-4/skovaka1@jhu.edu/code/UNCALLED/src/uncalled/models/r94_5mers.txt";
//const std::string PROBFN = "/home-4/skovaka1@jhu.edu/code/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

int main(int argc, char** argv) {

    std::string index(argv[1]), reads_fname(argv[2]);
    u32 nthreads = atoi(argv[3]);//, read_count = atoi(argv[4]);
    
    Params::init_map(index, MODEL,
                        22,    //seed_len
                        25,    //min_aln_len
                        0,     //min_rep_len
                        50,    //max_rep_copy
                        8,     //max_consec_stay
                        10000, //max_paths
                        30000, //max_events_proc
                        3,     //e/vt_winlen1
                        6,     //evt_winlen2
                        nthreads,     //nthreads
                        512,   //num_channels
                        1.4,   //evt_thresh1
                        9.0,   //evt_thresh2
                        0.2,   //evt_peak_height
                        30,    //evt_min_mean
                        150,   //evt_max_mean
                        0.5,   //max_stay_frac
                        -3.75, //min_seed_prob
                        7.00,  //min_mean_conf
                        2.25  //min_top_conf
                        );

    Timer t;

    Fast5Pool pool(reads_fname);
    u64 MAX_SLEEP = 100;

    while (!pool.all_finished()) {
        u64 t0 = t.get();
        for (Paf p : pool.update()) {
            p.print_paf();
        }
        u64 dt = t.get() - t0;
        if (dt < MAX_SLEEP) usleep(1000*(MAX_SLEEP - dt));
    }
}
