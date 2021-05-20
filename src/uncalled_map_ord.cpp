#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include "map_pool_ord.hpp"

void load_conf(int argc, char** argv, Conf &conf);

int main(int argc, char** argv) {
    std::cerr << "Loading conf\n";

    Conf conf(Conf::Mode::MAP_ORD);
    conf.set_mode_map_ord();

    load_conf(argc, argv, conf);

    MapPoolOrd pool(conf);

    Timer t;

    std::cerr << "Loading fast5s\n";

    pool.load_fast5s();

    u64 MAX_SLEEP = 100;

    std::cerr << t.lap() << " loaded\n";

    std::cerr << "Mapping\n";


    while (pool.running()) {
        u64 t0 = t.get();
        for (auto p : pool.update()) {
            std::get<2>(p).print_paf();
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

void load_conf(int argc, char** argv, Conf &conf) {
    int opt;
    std::string flagstr = "C:t:n:c:l:s:w:p:Sr";

    #ifdef DEBUG_OUT
    flagstr += "D:";
    #endif

    //parse flags
    while((opt = getopt(argc, argv, flagstr.c_str())) != -1) {
        switch(opt) {  

            FLAG_TO_CONF('t', atoi, threads)
            FLAG_TO_CONF('n', atoi, fast5_reader.max_reads)
            FLAG_TO_CONF('c', atoi, read_buffer.max_chunks)
            //FLAG_TO_CONF('l', std::string, fast5_reader.read_filter)
            FLAG_TO_CONF('p', std::string, mapper.idx_preset)

            #ifdef DEBUG_OUT
            FLAG_TO_CONF('D', std::string, mapper.dbg_prefix);
            #endif

            case 'l':
            conf.fast5_reader.load_read_filter(std::string(optarg));
            break;
   
            case 'C':
            std::cerr << "Conf: " << optarg << "\n";
            conf.load_toml(std::string(optarg));
            break;

            case 'r':
            std::cerr << "Setting RNA parameters\n";
            conf.set_r94_rna();
            break;

            case 'S':
            conf.fast5_reader.load_bc = true;
            conf.read_buffer.skip_notempl = true;
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

    POSITIONAL_TO_CONF(std::string, mapper.bwa_prefix)
    //POSITIONAL_TO_CONF(std::string, fast5_reader.fast5_list)
    if (i < argc) { 
        conf.fast5_reader.load_fast5_files(std::string(argv[i])); 
        i++; 
    } else { 
        std::cerr << "Error: must specify fast5_list\n"; \
        abort(); 
    } 
}
