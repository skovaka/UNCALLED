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
      fast5s_(conf.fast5_prms),
      scan_start_(0),
      is_running_(false),
      in_scan_(false) {

    float sample_rate = conf.get_sample_rate();
    time_coef_  = sample_rate / 1000;
    ej_time_    = PRMS.ej_time   * sample_rate;
    scan_time_  = PRMS.scan_time * sample_rate;
    max_chunks_ = conf.get_max_chunks();

    channels_.reserve(conf.num_channels);
    for (u32 c = 1; c <= conf.num_channels; c++) {
        channels_.emplace_back(c);
    }
}

bool ClientSim::load_from_files(const std::string &prefix) {
    std::string itvs_file = prefix + "_itvs.txt",
                gaps_file = prefix + "_gaps.txt",
                delays_file = prefix + "_delays.txt",
                reads_file = prefix + "_reads.txt";

    std::cerr << "Loading " << itvs_file << "\n";
    if (!load_itvs(itvs_file)) return false;

    std::cerr << "Loading " << gaps_file << "\n";
    if (!load_gaps(gaps_file)) return false;

    std::cerr << "Loading " << delays_file << "\n";
    if (!load_delays(delays_file)) return false;

    Timer t;

    std::cerr << "Loading reads\n";
    if (!load_reads(reads_file)) return false;
    std::cerr << "Loaded " << (t.get()/1000) << "\n";

    return true;
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

    infile >> ch >> i >> st >> en;
    while (!infile.eof()) {
        add_intv(ch, i, st, en);
        infile >> ch >> i >> st >> en;
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
    u32 len;

    infile >> ch >> i >> len;
    while (!infile.eof()) {
        add_gap(ch, i, len);
        infile >> ch >> i >> len;
    }

    return true;
}

bool ClientSim::load_delays(const std::string &fname) {
    if (fname.empty()) {
        std::cerr << "No delay file provided, using read delays\n";
        return false;
    }
    
    std::ifstream infile(fname);

    if (!infile.is_open()) {
        std::cerr << "Error: failed to open delay file\n";
        return false;
    }

    u16 ch, i;
    u32 len;

    infile >> ch >> i >> len;
    while (!infile.eof()) {
        infile >> ch >> i >> len;
        add_delay(ch, i, len);
    }

    return true;
}

void ClientSim::add_intv(u16 ch, u16 i, u32 st, u32 en) {
    channels_[ch-1].set_active(i, st, en);
}

void ClientSim::add_gap(u16 ch, u16 i, u32 len) {
    channels_[ch-1].add_gap(i, len);
}

void ClientSim::add_delay(u16 ch, u16 i, u32 len) {
    channels_[ch-1].add_delay(i, len);
}

void ClientSim::add_read(u16 ch, const std::string &id, u32 offs) {
    read_locs[id] = {ch, channels_[ch-1].reserve_read(), offs};
    fast5s_.add_read(id);
}

void ClientSim::add_fast5(const std::string &fname) {
    fast5s_.add_fast5(fname);
}


void ClientSim::load_fast5s() {
    u32 n = 0;
    while(!fast5s_.empty()) {
        ReadBuffer read = fast5s_.pop_read();
        ReadLoc r = read_locs[read.get_id()];

        read.set_channel(r.ch);
        channels_[r.ch-1].load_read(r.i, r.offs, read, max_chunks_);

        if (n % 1000 == 0) {
            std::cerr << n << " loaded\n";
        }
        n++;
    }
}


bool ClientSim::load_reads(const std::string &fname) {
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
    u32 offs;

    infile >> ch >> rd >> offs;
    while (!infile.eof()) {
        add_read(ch, rd, offs);
        infile >> ch >> rd >> offs;
    }


    return true;
}


bool ClientSim::run() {
    is_running_ = true;
    in_scan_ = false;
    timer_.reset();
    for (SimChannel &ch : channels_) {
        ch.start(0);
    }
    return true;
}

std::vector< std::pair<u16, Chunk> > ClientSim::get_read_chunks() {
    std::vector< std::pair<u16, Chunk> > ret; //TODO rename chunks?

    if (!is_running_) {
        return ret;
    }

    u64 time = get_time();

    //guilty till proven innocent
    bool intvs_ended = true; 
    bool next_intv = false;

    if (in_scan_) {
        if (time-scan_start_ >= scan_time_) {
            intvs_ended = in_scan_ = false;
            next_intv = true;
            std::cerr << time << " ending mux scan\n";
        } else {
            return ret;
        }
    }

    is_running_ = false;

    for (u16 c = 0; c < channels_.size(); c++) {
        SimChannel &ch = channels_[c];

        if (ch.is_dead()) continue;

        if (next_intv) {
            ch.next_intv(time);
            if (ch.is_dead()) continue;
        }

        is_running_ = true; 

        if (!ch.is_active(time)) {
            intvs_ended = ch.intv_ended(time) && intvs_ended;
            continue;
        }

        intvs_ended = false;

        while (ch.chunk_ready(time)) {
            ret.push_back( std::pair<u16, Chunk>(c+1, ch.next_chunk(time)));
        }
    }

    if (intvs_ended && !in_scan_) {
        std::cerr << time << " starting mux scan\n";
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

u32 ClientSim::unblock_read(u16 ch, u32 number) {
    if (get_number(ch) != number) { 
        return 0;
    }
    return channels_[ch-1].unblock(get_time(), ej_time_);
}

float ClientSim::get_time() {
    return (timer_.get() * time_coef_);
}

float ClientSim::get_runtime() {
    return timer_.get() / 1000;
}

bool ClientSim::is_running() {
    return is_running_;
}
