#include <iostream>
#include <unistd.h>
#include "conf.hpp"
#include "client_sim.hpp"
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


    ClientSim sim(conf);

    std::cerr << "Loading mappers\n";
    RealtimePool pool(conf);

    const u64 MAX_SLEEP = 100;

    std::cerr << "Starting simulation\n";
    sim.start();
    Timer t;

    std::vector<float> chunk_times(conf.get_num_channels(), t.get());

    while (sim.is_running()) {
        u64 t0 = t.get();

        u16 channel;
        u32 number;
        Paf paf;
        for (MapResult m : pool.update()) {
            std::tie(channel, number, paf) = m;
            float map_time = (t.get() - chunk_times[channel-1])/1000;
            if (paf.is_mapped()) {
                paf.set_float(Paf::Tag::EJECT, map_time); //TODO: match with realtime
                sim.unblock(channel, number);
                //aln.set_unblocked();TODO replace this
            } else {
                sim.stop_receiving_read(channel, number);
                paf.set_float(Paf::Tag::KEEP, map_time);
            }
            paf.print_paf();
            std::cout.flush();
        }

        std::vector<Chunk> chunks = sim.get_read_chunks();
        for (Chunk &ch : chunks) {
            if (pool.add_chunk(ch)) {
                chunk_times[ch.get_channel_idx()] = t.get();
            } else {
                std::cerr << "Error: failed to add chunk from " << ch.get_id() << std::endl;
            }
        }

        u64 dt = t.get() - t0;
        if (dt < MAX_SLEEP) usleep(1000*(MAX_SLEEP - dt));
    }
    std::cerr << "Reads aligned (" << (t.get() / 1000) << " sec)\n";

    pool.stop_all();
    std::cerr << "Done " << (t.get() / 1000) << "\n";
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
    while((opt = getopt(argc, argv, ":t:l:s:e:g:d:")) != -1) {
        switch(opt) {  

            FLAG_TO_CONF('l', std::string, read_list)
            FLAG_TO_CONF('g', std::string, gap_file)
            FLAG_TO_CONF('t', atoi, threads)
            FLAG_TO_CONF('s', atof, sim_start)
            FLAG_TO_CONF('e', atof, sim_end)
            FLAG_TO_CONF('d', atof, ej_delay)

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
