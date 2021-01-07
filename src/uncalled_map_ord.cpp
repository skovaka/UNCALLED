#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include "map_pool_ord.hpp"

void load_conf(int argc, char** argv, Conf &conf);

int main(int argc, char** argv) {
    std::cerr << "Loading conf\n";

    Conf conf(Conf::Mode::MAP_ORD);

    load_conf(argc, argv, conf);

    MapPoolOrd pool(conf);

    Timer t;

    pool.load_fast5s();

    u64 MAX_SLEEP = 100;

    std::cerr << t.lap() << " loaded\n";

    std::cerr << "Mapping\n";


    while (pool.running()) {
        u64 t0 = t.get();
        for (Paf p : pool.update()) {
            p.print_paf();
        }
        u64 dt = t.get() - t0;
        if (dt < MAX_SLEEP) usleep(1000*(MAX_SLEEP - dt));
    }

    std::cerr << t.lap() << " mapped\n";

    std::cerr << "Finishing\n";

    pool.stop();

}

#define FLAG_TO_CONF(C, T, F) { \
    case C: \
        conf.set_##F(T(optarg)); \
        break; \
}

#define POSITIONAL_TO_CONF(T, F) {\
    if (i < argc) { \
        conf.set_##F(T(argv[i])); \
        i++; \
    } else { \
        std::cerr << "Error: must specify flag " << #F << "\n"; \
        abort(); \
    } \
}

void load_conf(int argc, char** argv, Conf &conf) {
    int opt;
    std::string flagstr = "C:t:n:r:R:c:l:s:w:p:";

    #ifdef DEBUG_OUT
    flagstr += "D:";
    #endif

    //parse flags
    while((opt = getopt(argc, argv, flagstr.c_str())) != -1) {
        switch(opt) {  

            FLAG_TO_CONF('t', atoi, threads)
            FLAG_TO_CONF('n', atoi, max_reads)
            FLAG_TO_CONF('r', atoi, min_active_reads)
            FLAG_TO_CONF('R', atoi, max_active_reads)
            FLAG_TO_CONF('c', atoi, max_chunks)
            FLAG_TO_CONF('l', std::string, read_list)
            FLAG_TO_CONF('s', atof, win_stdv_min)
            FLAG_TO_CONF('w', atof, win_len)
            FLAG_TO_CONF('p', std::string, idx_preset)
            #ifdef DEBUG_OUT
            FLAG_TO_CONF('D', std::string, dbg_prefix);
            #endif

            case 'C':
            std::cerr << "Conf: " << optarg << "\n";
            conf.load_toml(std::string(optarg));
            break;

            case ':':  
            std::cerr << "Error: failed to load flag value\n";
            abort();

            case '?':  
            std::cerr << "Error: unknown flag\n";  
            abort();
        }
    }


    //parse positionals
    int i = optind;

    POSITIONAL_TO_CONF(std::string, bwa_prefix)
    POSITIONAL_TO_CONF(std::string, fast5_list)
}
