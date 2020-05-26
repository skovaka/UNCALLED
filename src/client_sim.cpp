/* MIT License
 *
 * Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "client_sim.hpp"

ClientSim::ClientSim(Conf &conf) :
      PRMS(conf.sim_prms),
      scan_start_(0),
      is_running_(false),
      in_scan_(false),
      channels_(conf.num_channels) {

    float sample_rate = conf.get_sample_rate();
    time_coef_  = PRMS.sim_speed * sample_rate / 1000;
    //start_samp_ = PRMS.sim_start * sample_rate;
    ej_delay_   = PRMS.ej_delay  * sample_rate;
    ej_time_    = PRMS.ej_time   * sample_rate;
    scan_time_  = PRMS.scan_time * sample_rate;

    //if (PRMS.sim_end > 0) end_samp_ = PRMS.sim_end * sample_rate;
    //else end_samp_ = INT_MAX;

    std::string prefix = conf.get_pat_prefix(),
                itvs_file = prefix + "_itvs.txt",
                gaps_file = prefix + "_gaps.txt",
               reads_file = prefix + "_reads.txt";

    std::cerr << "Loading intervals\n";
    load_itvs(itvs_file);

    std::cerr << "Loading gaps\n";
    load_gaps(gaps_file);

    std::cerr << "Loading reads\n";
    Fast5Reader fast5s(conf.fast5_prms);
    load_reads(reads_file, fast5s, Mapper::PRMS.max_chunks);
    std::cerr << "Loaded\n";
}


bool ClientSim::load_itvs(const std::string &fname) {
    if (fname.empty()) {
        std::cerr << "No gap file provided, using read gaps\n";
        return false;
    }
    
    std::ifstream infile(fname);

    if (!infile.is_open()) {
        std::cerr << "Error: failed to open gap file\n";
        return false;
    }

    u16 ch, i;
    u32 st, en;

    while (!infile.eof()) {
        infile >> ch >> i >> st >> en;
        channels_[ch-1].set_active(i, st, en);
    }

    return true;
}

bool ClientSim::load_gaps(const std::string &fname) {
    if (fname.empty()) {
        std::cerr << "No gap file provided, using read gaps\n";
        return false;
    }
    
    std::ifstream infile(fname);

    if (!infile.is_open()) {
        std::cerr << "Error: failed to open gap file\n";
        return false;
    }

    u16 ch, i;
    u32 ln;

    while (!infile.eof()) {
        infile >> ch >> i >> ln;
        channels_[ch-1].add_gap(i, ln);
    }

    return true;
}

bool ClientSim::load_reads(const std::string &fname, Fast5Reader &fast5s, u32 max_chunks) {
    std::unordered_map< std::string, std::pair<u16,u32> >  read_locs;
    std::vector<u32> ch_counts(channels_.size(), 0);

    if (fname.empty()) {
        std::cerr << "No read file provided\n";
        return false;
    }
    
    std::ifstream infile(fname);

    if (!infile.is_open()) {
        std::cerr << "Error: failed to open read file\n";
        return false;
    }

    u16 ch;
    std::string rd;

    while (!infile.eof()) {
        infile >> ch >> rd;
        fast5s.add_read(rd);
        read_locs[rd] = std::pair<u16,u32>(ch, ch_counts[ch-1]++);
    }

    for (u16 c = 0; c < ch_counts.size(); c++) {
        channels_[c].set_read_count(ch_counts[c]);
    }

    u32 n = 0;

    while(!fast5s.empty()) {
        ReadBuffer read = fast5s.pop_read();

        auto ch_i = read_locs[read.get_id()];
        u16 ch = ch_i.first;
        u32 i = ch_i.second;

        //std::cout << ch << " "
        //          << i << " " 
        //          << ch_counts[ch-1] << "\n";
        //std::cout.flush();

        read.set_channel(ch);
        channels_[ch-1].load_read(i, read, max_chunks);

        if (n % 1000 == 0) {
            std::cerr << n << " loaded\n";
        }
        n++;
    }

    return true;
}

void ClientSim::start() {
    is_running_ = true;
    in_scan_ = false;
    timer_.reset();
    u32 c = 0;
    for (SimChannel &ch : channels_) {
        if (ch.is_active(0)){ 
            //std::cerr << "INIT "
            //          << c << " "
            //          << ch.intv().gaps_.size() << "\n";
            ch.start(0);
            c++;
        }
    }
}



std::vector<Chunk> ClientSim::get_read_chunks() {
    std::vector<Chunk> ret; //TODO rename chunks?

    if (!is_running_) {
        return ret;
    }

    u64 time = get_time();

    //guilty till proven innocent
    is_running_ = false;
    bool intvs_ended = true; 
    bool next_intv = false;

    if (in_scan_) {
        if (time-scan_start_ >= scan_time_) {
            intvs_ended = in_scan_ = false;
            next_intv = true;
            std::cerr << "Ending mux scan\n";
        } else {
            return ret;
        }
    }


    for (u32 c = 0; c < channels_.size(); c++) {
        SimChannel &ch = channels_[c];

        if (ch.is_dead()) continue;

        if (next_intv) {
            ch.next_intv(time);
            if (ch.is_dead()) continue;
        }

        is_running_ = true; 

        //std::cerr << "Checking " << c << "\n";
        if (!ch.is_active(time)) {
            intvs_ended = ch.intv_ended(time) && intvs_ended;
            continue;
        }

        //std::cerr << c << " not dead\n";

        intvs_ended = false;

        while (ch.chunk_ready(time)) {
            //std::cerr << c << " kachunk\n";
            ret.push_back(ch.next_chunk(time));
        }
    }

    if (intvs_ended && !in_scan_) {
        std::cerr << "Starting mux scan\n";
        scan_start_ = time;
    }
    in_scan_ = intvs_ended;

    return ret;
}

u32 ClientSim::get_number(u16 ch) {
    return channels_[ch-1].read_number();
    //u32 c = channel-1;
    //if (channels_[c].empty()) return 0; 
    //auto &chunks = channels_[c].front().chunks;
    //if (chunks.empty()) return 0;
    //return chunks.front().get_number();
}

void ClientSim::stop_receiving_read(u16 ch, u32 number) {
    if (get_number(ch) == number) {
        channels_[ch-1].stop_receiving_read();
    }
}

void ClientSim::unblock(u16 ch, u32 number) {
    if (get_number(ch) != number) return;
    channels_[ch-1].unblock(get_time(), ej_delay_, ej_time_);
}

float ClientSim::get_time() {
    return (timer_.get() * time_coef_);
}


bool ClientSim::is_running() {
    return is_running_;
}
