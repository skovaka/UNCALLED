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

#ifndef FAST5_READER_HPP
#define FAST5_READER_HPP

#include <deque>
#include "fast5.hpp"
#include "fast5/hdf5_tools.hpp"
#include "util.hpp"
#include "timer.hpp"
#include "chunk.hpp"
#include "uncalled_opts.hpp"
#include "read_buffer.hpp" //TODO: separate ReadLoc so I don't have to do this


class Simulator {
    public:
    Simulator(const UncalledOpts &opts);
    void add_fast5s(const std::string &fname, u32 max_loaded);
    void start();

    std::vector<Chunk> get_read_chunks();
    void stop_receiving_read(u16 channel, u32 number);
    void unblock(u16 channel, u32 number);
    float get_time(u16 channel);
    bool is_running();


    private:
    u32 num_loaded_;
    u16 chunk_len_;
    float speed_;
    Timer timer_;
    bool is_running_;
    u64 tshift_;
    std::vector< u64 > chshifts_;
    std::vector< std::deque<Chunk> > chunks_;
};

#endif
