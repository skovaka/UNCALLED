#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include "map_pool.hpp"

bool load_conf(int argc, char** argv, Conf &conf);

int main(int argc, char** argv) {
    std::cerr << "Loading conf\n";

    Conf conf;

    if (!load_conf(argc, argv, conf)) {
        return 1;
    }

    MapPool pool(conf);

    u64 MAX_SLEEP = 100;

    std::cerr << "Mapping\n";

    Timer t;

    while (pool.running()) {
        u64 t0 = t.get();
        for (Paf p : pool.update()) {
            p.print_paf();
        }
        u64 dt = t.get() - t0;
        if (dt < MAX_SLEEP) usleep(1000*(MAX_SLEEP - dt));
    }

    std::cerr << "Finishing\n";

    pool.stop();

}

#define FLAG_TO_CONF(C, T, F) { \
    case C: \
        conf.F = T(optarg); \
        break; \
}

#define POSITIONAL_TO_CONF(T, F) {\
    if (i < argc) { \
        conf.F = T(argv[i]); \
        i++; \
    } else { \
        std::cerr << "Error: must specify flag " << #F << "\n"; \
        abort(); \
    } \
}

bool load_conf(int argc, char** argv, Conf &conf) {
    int opt;
    std::string flagstr = ":t:n:l:";

    #ifdef DEBUG_OUT
    flagstr += "D:";
    #endif

    while((opt = getopt(argc, argv, flagstr.c_str())) != -1) {

        switch(opt) {  

            FLAG_TO_CONF('t', atoi, threads)
            FLAG_TO_CONF('n', atoi, fast5_reader.max_reads)
            //FLAG_TO_CONF('l', std::string, fast5_reader.read_filter)

            case 'l':
            conf.fast5_reader.load_read_filter(std::string(optarg));
            break;

            #ifdef DEBUG_OUT
            FLAG_TO_CONF('D', std::string, mapper.dbg_prefix);
            #endif

            case ':':  
            std::cerr << "Error: failed to load flag value\n";  
            return false;

            case '?':  
            std::cerr << "Error: unknown flag\n";  
            return false;
        }
    }


    //parse positionals
    int i = optind;

    POSITIONAL_TO_CONF(std::string, mapper.bwa_prefix)
    if (i < argc) { 
        conf.fast5_reader.load_fast5_list(std::string(argv[i])); 
        i++; 
    } else { 
        std::cerr << "Error: must specify fast5_list\n"; \
        abort(); 
    } 

    return true;
}
