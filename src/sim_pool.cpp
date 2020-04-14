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

#include "sim_pool.hpp"

SimPool::SimPool(Conf &conf) :
      is_running_(false),
      tshift_(-1),
      chshifts_(conf.num_channels, 0),
      chunks_(conf.num_channels) {

    u32 max_chunks = Mapper::PRMS.max_chunks;

    std::cerr << "Loading read chunks\n";

    Fast5Reader fast5s(conf.fast5_prms);

    Timer t;

    u32 n = 0;

    while (!fast5s.empty()) {
        fast5s.fill_buffer();
        while (fast5s.buffer_size() > 0) {
            ReadBuffer r = fast5s.pop_read();
            r.get_chunks(chunks_[r.get_channel_idx()], max_chunks);
            n += 1;
        }
        std::cout << (u32) t.get() << "\t" << n
                  << "\tbuffered\n";
        std::cout.flush();
        t.reset();
    }

    std::cerr << "Sorting chunks on...\n";

    for (u32 i = 0; i < chunks_.size(); i++) {
        std::cerr << "\tchannel " << (i+1) << "\n";
        std::sort(chunks_[i].begin(), chunks_[i].end());
    }
    std::cout << (u32) t.get() << "\t" << n
              << "\tsorted\n";
    
}

//std::vector<Chunk> Simulator::get_read_chunks() {
//    std::vector<Chunk> ret;
//
//    u64 time = (timer_.get() * speed_) + tshift_;
//
//    if (time > PRMS.end*PARAMS.sample_rate) {
//        is_running_ = false;
//        return ret;
//    }
//    
//    is_running_ = false;
//
//    for (u16 c = 0; c < chunks_.size(); c++) {
//        if (chunks_[c].empty()) continue;
//        is_running_ = true; 
//
//        //Find first chunk that ends after current time
//        //Ideally will be second chunk, unless we missed some
//        u16 i = 0;
//        for (; i < chunks_[c].size() && 
//               chunks_[c][i].get_start()+chunk_len_ < time+chshifts_[c]; i++);
//
//        //Skip if first chunk, otherwise add previous chunk
//        if (i-- == 0) {
//            continue;
//        }
//
//        chunks_[c][i].set_start(time-chunk_len_);
//        ret.push_back(chunks_[c][i]);
//        //TODO: try code below
//        //ret.push_back(Chunk());
//        //ret.back().swap(chunks_[c][i]);
//
//        //Remove all finished chunks
//        for (; i < chunks_[c].size(); i--) {
//            chunks_[c].pop_front();
//        }
//
//        if (chunks_[c].empty()) {
//            std::cout << "# " << c << " empty\n";
//        }
//    }
//
//    return ret;
//}
//
//void Simulator::stop_receiving_read(u16 channel, u32 number) {
//    channel--;
//    while (!chunks_[channel].empty() && 
//            chunks_[channel][0].get_number() == number) {
//        chunks_[channel].pop_front();
//    }
//}
//
//void Simulator::unblock(u16 channel, u32 number) {
//    channel--;
//    u64 t0 = chunks_[channel].front().get_start();
//    u64 end = t0;
//    while (!chunks_[channel].empty() && 
//            chunks_[channel][0].get_number() == number) {
//        end = chunks_[channel].front().get_end();
//        chunks_[channel].pop_front();
//    }
//    chshifts_[channel] += end - t0 - (PARAMS.sim_gaps*PARAMS.sample_rate);
//}
//
//float Simulator::get_time(u16 channel) {
//    channel--;
//    return (timer_.get() * speed_) + tshift_ + chshifts_[channel];
//}
//
//void Simulator::start() {
//    speed_ = PARAMS.sample_rate * PARAMS.sim_speed / 1000.0;
//    is_running_ = true;
//    timer_.reset();
//}
//
//bool Simulator::is_running() {
//    return is_running_;
//}
