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
    sim.run();
    Timer t;

    std::vector<float> chunk_times(conf.get_num_channels(), t.get());
    std::vector<u32> unblocked(conf.get_num_channels(), 0);

    bool deplete = conf.get_realtime_mode() == RealtimeParams::Mode::DEPLETE;

    std::cerr << "Starting " << deplete << "\n";

    while (sim.is_running()) {
        u64 t0 = t.get();

        u16 channel;
        u32 number;
        Paf paf;
        for (MapResult m : pool.update()) {
            std::tie(channel, number, paf) = m;
            float map_time = (t.get() - chunk_times[channel-1])/1000;

            if (paf.is_ended()) {
                paf.set_float(Paf::Tag::ENDED, map_time);
                sim.stop_receiving_read(channel, number);

            } else if ((paf.is_mapped() && deplete) || (!paf.is_mapped() && !deplete)) {

                u32 delay = sim.unblock(channel, number);
                paf.set_float(Paf::Tag::EJECT, map_time); 
                paf.set_int(Paf::Tag::DELAY, delay); 

                unblocked[channel-1] = number;

            } else {
                sim.stop_receiving_read(channel, number);
                paf.set_float(Paf::Tag::KEEP, map_time);
            }
            paf.print_paf();
            std::cout.flush();
        }

        for (auto &r : sim.get_read_chunks()) {
            Chunk &ch = r.second;
            if (unblocked[ch.get_channel_idx()] == ch.get_number()) {
                std::cout << "# recieved chunk from " 
                          << ch.get_id() 
                          << " after unblocking\n";
                continue;
            } else if (pool.add_chunk(ch)) {
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
    while((opt = getopt(argc, argv, ":t:l:s:g:c:p:de")) != -1) {
        switch(opt) {  

            FLAG_TO_CONF('l', std::string, read_list)
            FLAG_TO_CONF('s', std::string, sim_prefix)
            FLAG_TO_CONF('p', std::string, index_preset)
            FLAG_TO_CONF('t', atoi, threads)
            FLAG_TO_CONF('c', atoi, max_chunks)

            case 'd':
                conf.set_realtime_mode(RealtimeParams::Mode::DEPLETE);
                break;

            case 'e':
                conf.set_realtime_mode(RealtimeParams::Mode::ENRICH);
                break;

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
