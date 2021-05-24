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

#include <thread>
#include <chrono>
#include "pdqsort.h"
#include "map_pool_ord.hpp"

MapPoolOrd::MapPoolOrd(Config &config) : 
        RealtimePool(config),
        fast5s_(config.fast5_reader),
        channels_empty_(false) {

    config.set_mode_map_ord();
    config.export_static();

    channels_.resize(config.get_num_channels());
    chunk_idx_.resize(config.get_num_channels());

    #ifdef PYDEBUG
    meta_hold_.resize(config.get_num_channels(), 0);
    #endif
}

void MapPoolOrd::add_fast5(const std::string &fname) {
    fast5s_.add_fast5(fname);
}

void MapPoolOrd::add_read(const std::string &id) {
    fast5s_.add_read(id);
}

void MapPoolOrd::load_fast5s() {
    while(!fast5s_.empty()) {
        //ReadBuffer read = static_cast<ReadBuffer>(fast5s_.next_read());
        auto read = fast5s_.next_read();
        channels_[read.get_channel_idx()].push_back(read);
    }

    for (auto &ch : channels_) {
        pdqsort(ch.begin(), ch.end());
    }
}

std::vector<MapResult> MapPoolOrd::update() {

    channels_empty_ = true;

    for (u32 i = 0; i < channels_.size(); i++) {
        if (channels_[i].empty()) continue;
        channels_empty_ = false;

        #ifdef PYDEBUG
        if (meta_hold_[i]) continue;
        #endif

        //Skip the read if the mapper has finished
        //It will be removed from the queue below
        if (is_read_finished(channels_[i].front())) {
            continue;
        }

        //Get next chunk
        auto &r = channels_[i].front();
        
        //TODO turn chunk into slice of ReadBuffer/Fast5Read
        auto chunk = r.get_chunk(chunk_idx_[i]); 

        //Try adding to pool
        //If sucessfful, move to next chunk
        if (try_add_chunk(chunk)) {
            chunk_idx_[i]++;
        }
    }

    auto ret = RealtimePool::update();

    //Get mapping results
    for (auto &m : ret) {
        u16 channel = std::get<0>(m);
        #ifndef PYDEBUG
        end_read(channel);
        #else
        meta_hold_[channel-1] = true;
        #endif
    }

    if (active_count() < PRMS.min_active_reads) {
        stop_all();
        for (auto &chs : channels_) chs.clear();
        channels_empty_ = true;
    }

    return ret;
}

void MapPoolOrd::end_read(u16 channel) {
    channels_[channel-1].pop_front();
    chunk_idx_[channel-1] = 0;
    #ifdef PYDEBUG
    meta_hold_[channel-1] = false;
    #endif
}

Fast5Read &MapPoolOrd::get_read(u16 channel) {
    return channels_[channel-1].front();
}

bool MapPoolOrd::is_read_finished(const ReadBuffer &r) {
    auto &mapper = get_mapper(r.get_channel());
    return (mapper.finished() && 
            mapper.get_read().get_number() == r.get_number());
}

bool MapPoolOrd::try_add_chunk(ReadBuffer &chunk) {
    u16 ch = chunk.get_channel_idx();
    auto &mapper = get_mapper(ch+1);

    //Chunk is empty if all read chunks were output
    if (chunk.empty()) {

        //TODO probably dont need
        #ifdef PYDEBUG
        if (meta_hold_[ch]) return false;
        #endif

        //Give up if previous chunk done mapping
        if (mapper.chunk_mapped() && !mapper.finished()) {
            mapper.request_reset();
        }
        return false;
    }

    //Start new read if mapper inactive
    if (mapper.get_state() == Mapper::State::INACTIVE) {
        mapper.new_read(chunk);
        active_queue_.push_back(ch);
        return true;

    } else if (mapper.get_read().number_ == chunk.get_number()) {

        //Don't add if previous chunk is still mapping
        if (!mapper.chunk_mapped()) {
            return false;
        }

        return mapper.add_chunk(chunk);
    }

    return false;
}

bool MapPoolOrd::running() {
    return !(channels_empty_ && all_finished());
}

void MapPoolOrd::stop() {
    return stop_all(); //TODO rename stop in realtime
}

