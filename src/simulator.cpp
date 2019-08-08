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

#include "simulator.hpp"
#include "params.hpp"

Simulator::Simulator() 
    : chunk_len_(PARAMS.chunk_len),
      speed_(PARAMS.sample_rate * PARAMS.sim_speed / 1000),
      is_running_(false),
      tshift_(-1),
      chshifts_(PARAMS.num_channels, 0),
      chunks_(PARAMS.num_channels) {}

void Simulator::add_fast5s(const std::string &fname, u32 max_loaded) {
    std::deque<ReadBuffer> reads;
    load_fast5s(fname, reads, max_loaded);
    std::cerr << "Loaded " << reads.size() << " reads\n";
    //std::sort(reads.begin(), reads.end());

    tshift_ = reads.front().start_sample_;

    for (const ReadBuffer &r : reads) {
        r.get_chunks(chunks_[r.get_channel_idx()], chunk_len_);
    }
}

std::vector<Chunk> Simulator::get_read_chunks() {
    std::vector<Chunk> ret;

    if (timer_.get() / 1000.0 > (PARAMS.sim_en - PARAMS.sim_st+1)) {
        is_running_ = false;
        return ret;
    }

    u64 time = (timer_.get() * speed_) + tshift_;
    
    is_running_ = false;

    for (u16 c = 0; c < chunks_.size(); c++) {
        if (chunks_[c].empty()) continue;
        is_running_ = true; 

        //Find first chunk that ends after current time
        //Ideally will be second chunk, unless we missed some
        u16 i = 0;
        for (; i < chunks_[c].size() && 
               chunks_[c][i].get_start()+chunk_len_ < time+chshifts_[c]; i++);

        //Skip if first chunk, otherwise add previous chunk
        if (i-- == 0) {
            continue; 
        }

        ret.push_back(chunks_[c][i]);
        //TODO: try code below
        //ret.push_back(Chunk());
        //ret.back().swap(chunks_[c][i]);

        //Remove all finished chunks
        for (; i < chunks_[c].size(); i--) {
            chunks_[c].pop_front();
        }
    }

    return ret;
}

void Simulator::stop_receiving_read(u16 channel, u32 number) {
    channel--;
    while (!chunks_[channel].empty() && 
            chunks_[channel][0].get_number() == number) {
        chunks_[channel].pop_front();
    }
}

void Simulator::unblock(u16 channel, u32 number) {
    channel--;
    u64 t0 = chunks_[channel].front().get_start();
    u64 end = t0;
    while (!chunks_[channel].empty() && 
            chunks_[channel][0].get_number() == number) {
        end = chunks_[channel].front().get_end();
        chunks_[channel].pop_front();
    }
    chshifts_[channel] += end - t0 - (0.1*PARAMS.sample_rate);
}

float Simulator::get_time(u16 channel) {
    channel--;
    return (timer_.get() * speed_) + tshift_ + chshifts_[channel];
}

void Simulator::start() {
    speed_ = PARAMS.sample_rate * PARAMS.sim_speed / 1000.0;
    is_running_ = true;
    timer_.reset();
}

bool Simulator::is_running() {
    return is_running_;
}
