#include <iostream>
#include <unistd.h>
#include "conf.hpp"
#include "sim_pool.hpp"
#include "realtime_pool.hpp"
#include "read_buffer.hpp"

const std::string CONF_DIR(std::getenv("UNCALLED_CONF")),
                  DEF_MODEL = CONF_DIR + "/r94_5mers.txt",
                  DEF_CONF = CONF_DIR + "/defaults.toml";

bool load_conf(int argc, char** argv, Conf &conf);

int main(int argc, char** argv) {
    std::cerr << "Loading conf\n";

    Conf conf(DEF_CONF);
    conf.set_kmer_model(DEF_MODEL);

    if (!load_conf(argc, argv, conf)) {
        return 1;
    }

    Timer t;

    SimPool sim(conf);
    //sim.add_fast5s(reads_fname, read_count);

    //std::cerr << "Chunks split (" << (t.lap() / 1000) << " sec), aligning\n";
    //std::cerr.flush();

    //RealtimePool pool;

    //const u64 MAX_SLEEP = 100;

    //sim.start();
    //while (sim.is_running() || !pool.all_finished()) {
    //    u64 t0 = t.get();
    //    
    //    //if (!chunks.empty()) {
    //    //    std::cout << t.get() << "\t" << chunks.size() << " chunks\n";
    //    //}

    //    u16 channel;
    //    u32 number;
    //    Paf paf;
    //    for (MapResult m : pool.update()) {
    //        std::tie(channel, number, paf) = m;
    //        if (paf.is_mapped()) {
    //            paf.set_float(Paf::Tag::UNBLOCK, sim.get_time(channel));
    //            sim.unblock(channel, number);
    //            //aln.set_unblocked();TODO replace this
    //        } else {
    //            sim.stop_receiving_read(channel, number);
    //            paf.set_float(Paf::Tag::KEEP, sim.get_time(channel));
    //        }
    //        paf.print_paf();
    //    }

    //    std::vector<Chunk> chunks = sim.get_read_chunks();
    //    for (Chunk &ch : chunks) {
    //        if (!pool.add_chunk(ch)) {
    //            std::cerr << "Error: failed to add chunk from " << ch.get_id() << std::endl;
    //        }
    //    }

    //    u64 dt = t.get() - t0;
    //    if (dt < MAX_SLEEP) usleep(1000*(MAX_SLEEP - dt));
    //}

    //pool.stop_all();
    //std::cerr << "Reads aligned (" << (t.lap() / 1000) << " sec)\n";
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
        std::cerr << "Error: must specify " << #F << "\n"; \
        return false; \
    } \
}

bool load_conf(int argc, char** argv, Conf &conf) {
    int opt;

    //parse flags
    while((opt = getopt(argc, argv, ":t:n:l:")) != -1) {
        switch(opt) {  

            FLAG_TO_CONF('t', atoi, threads)
            FLAG_TO_CONF('l', std::string, read_list)

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

    POSITIONAL_TO_CONF(std::string, bwa_prefix)
    POSITIONAL_TO_CONF(std::string, fast5_list)

    return true;
}
