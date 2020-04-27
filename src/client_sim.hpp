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

#ifndef SIM_POOL_HPP
#define SIM_POOL_HPP

#include <deque>
#include "fast5.hpp"
#include "fast5/hdf5_tools.hpp"
#include "util.hpp"
#include "chunk.hpp"
#include "read_buffer.hpp" 
#include "fast5_reader.hpp" 
#include "conf.hpp" 

class ClientSim {
    public:
    ClientSim(Conf &c);
    std::vector< std::deque<float> > load_gaps(const std::string &fname, float sim_st);
    void add_fast5s(const std::string &fname, u32 max_loaded);
    void start();

    std::vector<Chunk> get_read_chunks();
    void stop_receiving_read(u16 channel, u32 number);
    void unblock(u16 channel, u32 number);
    float get_time();
    u32 get_number(u16 channel);
    bool is_running();

    private:

    class SimRead {
        public:
        SimRead(const ReadBuffer &r, u32 max_chunks);
        SimRead(float terminal_gap);

        bool in_gap(u64 t);
        bool in_read(u64 t);
        bool finished(u64 t);
        bool chunk_ready(u64 t);
        void make_terminal();
        bool is_terminal();
        u64 get_end();
        void set_start(u64 t);

        std::deque<Chunk> chunks;
        u64 start, gap, duration;
    };

    friend bool operator< (const SimRead &r1, const SimRead &r2);

    SimParams PRMS;
    float time_coef_; //TODO: make const?
    u32 start_samp_, end_samp_, ej_time_, ej_delay_;

    u32 num_loaded_;
    bool is_running_;

    //u64 tshift_;
    //std::vector< u64 > chshifts_;

    Timer timer_;

    std::vector< std::deque<SimRead> > channels_;
};

bool operator< (const ClientSim::SimRead &r1, const ClientSim::SimRead &r2);
#endif
