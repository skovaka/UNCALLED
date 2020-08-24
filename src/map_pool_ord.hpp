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

#ifndef MAP_POOL_ORD_HPP
#define MAP_POOL_ORD_HPP

#include <thread>
#include <vector>
#include <deque>
#include <unordered_set>
#include "conf.hpp"
#include "realtime_pool.hpp"
#include "fast5_reader.hpp"

class MapPoolOrd {
    public:

    MapPoolOrd(Conf &conf);

    void add_fast5(const std::string &fname);
    void add_read(const std::string &id);
    void load_fast5s();

    std::vector<Paf> update();
    void stop();
    bool running();

    private:
    Fast5Reader fast5s_;
    RealtimePool pool_;

    u32 active_tgt_;

    using ChQueue = std::deque<ReadBuffer>;
    std::vector<ChQueue> channels_;
    std::vector<u32> chunk_idx_;

    bool channels_empty_;
};


#endif
