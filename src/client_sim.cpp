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
      is_running_(false),
      channels_(conf.num_channels) {

    float sample_rate = conf.get_sample_rate();
    time_coef_  = PRMS.sim_speed * sample_rate / 1000;
    start_samp_ = PRMS.sim_start * sample_rate;
    start_samp_ = PRMS.sim_start * sample_rate;
    ej_delay_   = PRMS.ej_delay  * sample_rate;
    ej_time_    = PRMS.ej_time   * sample_rate;

    if (PRMS.sim_end > 0) end_samp_ = PRMS.sim_end * sample_rate;
    else end_samp_ = INT_MAX;

    u32 max_chunks = Mapper::PRMS.max_chunks;

    auto gaps = load_gaps(conf.get_gap_file(), PRMS.sim_start);

    std::cerr << "Loading read chunks\n";

    Fast5Reader fast5s(conf.fast5_prms);

    Timer t;

    u32 n = 0;

    while (!fast5s.empty()) {
        fast5s.fill_buffer();
        while (fast5s.buffer_size() > 0) {
            ReadBuffer read = fast5s.pop_read();
            if (read.get_start() < start_samp_) continue;

            u32 c = read.get_channel_idx();
            channels_[c].emplace_back(read, max_chunks);

            n += 1;
        }

        std::cerr << (u32) t.get() << "\t" << n
                  << "\tbuffered\n";
        std::cerr.flush();
        t.reset();
    }

    std::cerr << "Sorting\n";

    for (u32 i = 0; i < channels_.size(); i++) {
        if (channels_[i].empty()) continue;

        std::sort(channels_[i].begin(), channels_[i].end());

        if (gaps[i].empty()) {
            u64 prev_end = 0;
            for (SimRead &r : channels_[i]) {
                r.gap = r.start - prev_end;
                prev_end = r.start + r.duration;
            }
            channels_[i].back().make_terminal();
            std::cerr << "No more gaps\n";

        } else {
            for (SimRead &r : channels_[i]) {
                //TODO: temporary 3.0
                r.gap = (gaps[i].front()) * sample_rate;
                gaps[i].pop_front();
                if (gaps[i].empty()) break;
            }

            if (!gaps[i].empty()) {
                channels_[i].emplace_back(gaps[i].front() * sample_rate); //add terminal gap
            } else {
                std::cerr << "PROBLEM " << (i+1) << "\n";
            }
        }

        channels_[i].front().set_start(0);
    }

    std::cerr << (u32) t.get() << "\t" << n
              << "\tsorted\n";
}

std::vector< std::deque<float> > ClientSim::load_gaps(const std::string &fname, float sim_st) {
    std::vector< std::deque<float> > ret(channels_.size());

    if (fname.empty()) {
        std::cerr << "No gap file provided, using read gaps\n";
        return ret;
    }
    
    std::ifstream infile(fname);

    if (!infile.is_open()) {
        std::cerr << "Error: failed to open gap file\n";
        return ret;
    }

    std::cerr << "Loading gap file\n";

    u32 ch, mx1, mx2;
    float st, ln;

    //TODO: filter mux changes?
    //note: assuming already sorted, ignoring start
    while (!infile.eof()) {
        infile >> ch >> st >> ln >> mx1 >> mx2;
        if (st + ln >= sim_st) {
            ret[ch-1].push_back(ln);
        }
    }

    return ret;
}

void ClientSim::start() {
    is_running_ = true;
    timer_.reset();
}


std::vector<Chunk> ClientSim::get_read_chunks() {
    std::vector<Chunk> ret;

    if (!is_running_) {
        return ret;
    }

    u64 time = get_time();

    if (time > end_samp_) {
        is_running_ = false;
        return ret;
    }
    
    is_running_ = false;

    u16 channel = 0;
    for (auto &reads : channels_) {
        channel++;

        if (reads.empty()) continue;
        is_running_ = true; 

        if (reads.front().finished(time)) {
            if (!reads.front().chunks.empty()) {
                std::cout << "# skipping " 
                          << reads.front().chunks.size() 
                          << " chunks\n";
            }

            if (reads.front().is_terminal()) {
                std::cout << "# simulation finished (channel " 
                          << channel << " empty) "
                          << reads.front().gap << " "
                          << reads.front().get_end()  << "\n"; 
                is_running_ = false;
                return ret;
            }

            u64 post_time = time - reads.front().get_end();
            reads.pop_front();

            if (reads.empty()) {
                std::cout << "# sim channel " 
                          << channel 
                          << " empty\n";
                continue;
            }

            reads.front().set_start(time - post_time);
        }

        SimRead &read = reads.front();

        if (read.in_read(time)) {
            Chunk c;
            while (read.chunk_ready(time)) {
                c = read.chunks.front();
                read.chunks.pop_front();
            }
            if (!c.empty()) ret.push_back(c);
        }

    }

    return ret;
}

u32 ClientSim::get_number(u16 channel) {
    u32 c = channel-1;
    if (channels_[c].empty()) return 0; 

    auto &chunks = channels_[c].front().chunks;
    if (chunks.empty()) return 0;

    return chunks.front().get_number();
}

void ClientSim::stop_receiving_read(u16 channel, u32 number) {
    if (get_number(channel) == number) {
        channels_[channel-1].front().chunks.clear();
    }
}

void ClientSim::unblock(u16 channel, u32 number) {
    if (get_number(channel) != number) return;

    auto r = channels_[channel-1].begin();

    u64 t = get_time();

    if (r->in_read(t)) {
        u32 new_dur = get_time() - (r->start + r->gap) + ej_delay_;
        if (new_dur < r->duration) r->duration = new_dur;
        r->chunks.clear();
    } else if (!r->finished(t)) {
        std::cerr << "I DON'T THINK THIS SHOULD HAPPEN?\n";
    }

    if (r != channels_[channel-1].end()) {
        r++;
        r->gap += ej_time_;
    }
    
    //TODO: increase next read gap by ej_samps_

    //auto &reads = channels_[channel-1];
    //reads.pop_front();

    //if (!reads.empty()) {
    //    reads.front().set_start(get_time() + ej_samps_);
    //}
}

float ClientSim::get_time() {
    return (timer_.get() * time_coef_);
}


bool ClientSim::is_running() {
    return is_running_;
}

ClientSim::SimRead::SimRead(const ReadBuffer &r, u32 max_chunks) :
    start (r.get_start()),
    duration (r.get_duration()) {
    r.get_chunks(chunks, max_chunks, false);
}

ClientSim::SimRead::SimRead(float terminal_gap) :
    gap (terminal_gap),
    duration (0) {}

bool ClientSim::SimRead::in_gap(u64 t) {
    return t >= start && t < start+gap;
}

bool ClientSim::SimRead::in_read(u64 t) {
    return t >= start+gap && t < get_end();
}

bool ClientSim::SimRead::finished(u64 t) {
    return t >= get_end();
}

u64 ClientSim::SimRead::get_end() {
    return start+gap+duration;
}

bool ClientSim::SimRead::chunk_ready(u64 t) {
    return !chunks.empty() && t >= chunks.front().get_end();
}

void ClientSim::SimRead::make_terminal() {
    duration = 0;
    chunks.clear();
}

bool ClientSim::SimRead::is_terminal() {
    return duration == 0;
}

void ClientSim::SimRead::set_start(u64 t) {
    start = t;
    u64 i = start+gap;
    for (auto &c : chunks) {
        c.set_start(i);
        i += c.size();
    }
}

bool operator< (const ClientSim::SimRead &r1, const ClientSim::SimRead &r2) {
    return r1.start < r2.start;
}
