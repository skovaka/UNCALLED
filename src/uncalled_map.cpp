#include <iostream>
#include <unistd.h>
#include "fast5_pool.hpp"

//const std::string MODEL =  "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers.txt";
//const std::string PROBFN = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

const std::string MODEL = "/home/skovaka/code/UNCALLED/src/uncalled/models/r94_5mers.txt";
const std::string PROBFN = "/home/skovaka/code/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

//const std::string MODEL = "/home-4/skovaka1@jhu.edu/code/UNCALLED/src/uncalled/models/r94_5mers.txt";
//const std::string PROBFN = "/home-4/skovaka1@jhu.edu/code/UNCALLED/src/uncalled/models/r94_5mers_threshs.txt";

int main(int argc, char** argv) {

    std::string index(argv[1]), reads_fname(argv[2]), conf_fname(argv[3]);

    Conf conf(conf_fname);
    
    Timer t;

    Fast5Pool pool(conf);
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
