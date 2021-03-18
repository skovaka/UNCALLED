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

    Conf conf;

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

                u32 delay = sim.unblock_read(channel, number);
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
            auto &chunk = r.second;
            if (unblocked[chunk.get_channel_idx()] == chunk.get_number()) {
                std::cout << "# recieved chunk from " 
                          << chunk.get_id() 
                          << " after unblocking\n";
                continue;
            } else if (pool.add_chunk(chunk)) {
                chunk_times[chunk.get_channel_idx()] = t.get();
            } else {
                std::cerr << "Error: failed to add chunk from " << chunk.get_id() << std::endl;
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
    std::string flagstr = ":t:c:p:de";

    #ifdef DEBUG_OUT
    flagstr += "D:";
    #endif

    //parse flags
    while((opt = getopt(argc, argv, flagstr.c_str())) != -1) {
        switch(opt) {  

            FLAG_TO_CONF('t', atoi, threads)
            FLAG_TO_CONF('p', std::string, mapper.idx_preset)
            FLAG_TO_CONF('c', atoi, read_buffer.max_chunks)

            #ifdef DEBUG_OUT
            FLAG_TO_CONF('D', std::string, mapper.dbg_prefix);
            #endif

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
