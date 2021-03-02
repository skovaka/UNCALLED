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

#ifndef _INCL_READ_BUFFER
#define _INCL_READ_BUFFER

/*TODO 
 * need stronger distinction between full read with chunks at offsets
 * vs rolling chunk input from live mode
 * vs simulated chunks (which we could avoid copying)
 * basically need to reconcile realtime, client_sim, and map_ord
*/

#include <string>
#include <vector>
#include <unordered_set>
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#endif

class ReadBuffer {
    public:

    typedef struct {
        u16 num_channels;
        float bp_per_sec;
        float sample_rate;
        float chunk_time;
        u32 start_chunk;
        u32 max_chunks;
        bool seq_fwd; //If true reads start at 5', end at 3' (DNA)
                      //If false start=3', end=5' (RNA)
        bool skip_notempl;

        float bp_per_samp() {
            return bp_per_sec / sample_rate;
        }

        u16 chunk_len() {
            return (u16) (chunk_time * sample_rate);
        }
    } Params;

    static Params PRMS;

    ReadBuffer();

    ReadBuffer(const std::string &id, u16 channel, u32 number, u64 start_time, 
          const std::string &dtype, const std::string &raw_str);
    ReadBuffer(const std::string &id, u16 channel, u32 number, u64 start_time, 
          const std::vector<float> &raw_data, u32 raw_st, u32 raw_len);

    //Returns read ID (name)
    std::string get_id() const;

    //Returns channel number and index 
    u16 get_channel() const;
    void set_channel(u16 ch); //{channel_idx_ = ch-1;}
    u16 get_channel_idx() const;

    u32 get_number() const; //{return number_;}

    //Returns start sample of read
    u64 get_start() const;

    //Sets the start sample
    void set_start(u64 start);

    //Returns last known sample of read
    u64 get_end() const;

    //Returns size of buffered signal in samples
    u64 size() const;
    
    //Returns total known size of read
    //if chunk_count()==1 
    u64 get_full_duration() const; 
    void set_full_duration(u64 raw_len_); 

    //Returns the signal buffer
    const std::vector<float> &get_signal() const;

    //Returns sample at specified index
    float &operator[](u32 i);

    void swap(ReadBuffer &r);

    //Clears the signal buffer
    void clear();
    
    //Returns true if no signal is in buffer
    bool empty() const;

    //Returns number of chunks loaded in buffer
    u32 get_chunk_count() const;

    //Returns true if the maximal number of chunks are in the buffer
    bool chunks_maxed() const ;

    //Returns the ith chunk in the buffer
    //only used if chunk_count() > 1
    ReadBuffer get_chunk(u32 i) const;

    //TODO make simulator just read from ReadBuffer
    //only used if chunk_count() > 1
    u32 get_chunks(std::vector<ReadBuffer> &chunk_queue, bool real_start=true, u32 offs=0) const;

    //Updates the read chunk
    //only used if chunk_count() == 1
    bool add_next_chunk(ReadBuffer &r);

    #ifdef PYBIND

    #define PY_READ_METH(P) c.def(#P, &ReadBuffer::P);
    #define PY_READ_RPROP(P) c.def_property_readonly(#P, &ReadBuffer::get_##P);
    #define PY_READ_PRM(P) p.def_readwrite(#P, &ReadBuffer::Params::P);

    static void pybind_defs(pybind11::class_<ReadBuffer> &c) {
        c.def(pybind11::init<ReadBuffer>());
        PY_READ_METH(empty);
        PY_READ_RPROP(id);
        PY_READ_RPROP(start);
        PY_READ_RPROP(end);
        PY_READ_RPROP(full_duration);
        PY_READ_RPROP(channel);
        PY_READ_RPROP(number);
        PY_READ_RPROP(signal);

        c.def("__len__", &ReadBuffer::size);
        c.def("__getitem__", &ReadBuffer::operator[]);

        pybind11::class_<Params> p(c, "Params");
        PY_READ_PRM(num_channels);
        PY_READ_PRM(bp_per_sec);
        PY_READ_PRM(sample_rate);
        PY_READ_PRM(chunk_time);
        PY_READ_PRM(start_chunk);
        PY_READ_PRM(max_chunks);
        PY_READ_PRM(load_bc);
        PY_READ_PRM(skip_notempl);
    }

    #endif

    //Source source_;
    std::string id_;
    u16 channel_idx_;
    u32 number_;
    u64 start_sample_, full_duration_; //TODO no raw len
    std::vector<float> signal_; //TODO store one
    u16 chunk_count_; //TODO derive from chunk_len?
    bool single_chunk_;

    friend bool operator< (const ReadBuffer &r1, const ReadBuffer &r2);
};

bool operator< (const ReadBuffer &r1, const ReadBuffer &r2);

#endif
