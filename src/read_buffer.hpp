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
 * Refactor into Paf and SigBuffer
 * Paf replaces python Paf
 *      seprate from bufs, maybe depends on them
 *      put cigar (+cs string) parsing here
 * SigBuffer stores signal, time, channel, read num, etc
 *      static link to sample rate, chunk len?
 *      merge with Chunk, don't duplicate data
 * 
 * Fast5Read : SigBuffer subclass?
 *      keep fast5 parsing stuff seperate
 *      optionally load bc data (moves)
*/

#include <string>
#include <vector>
#include <unordered_set>
//#include <fast5/hdf5_tools.hpp>
#include "hdf5_tools.hpp"
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
        u32 max_chunks;

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

    //TODO move to subclass, maybe friend to fast5reader?
    ReadBuffer(const hdf5_tools::File &file, const std::string &raw_path, const std::string &ch_path, const std::string &seg_path="");
    

    //Returns true if buffer contains a single chunk
    bool single_chunk() const;

    //Returns true if no signal is in buffer
    bool empty() const;

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
    //equal to size() if !single_chunk()
    u64 full_duration() const; 
    void set_full_duration(u64 raw_len_); //_raw_len

    //Returns the signal buffer
    const std::vector<float> &get_signal() const; //{return full_signal_;}

    //Returns sample at specified index
    float &operator[](u32 i);// {return full_signal_[i];}

    void swap(ReadBuffer &r);

    void clear();

    //Returns number of chunks loaded in buffer
    u32 chunk_count() const;

    //Returns true if the maximal number of chunks are in the buffer
    bool chunks_maxed() const ;

    //Returns the ith chunk in the buffer
    //only relavent if !single_chunk()
    ReadBuffer get_chunk(u32 i) const;

    //TODO make simulator just read from ReadBuffer
    //only relavent if !single_chunk()
    u32 get_chunks(std::vector<ReadBuffer> &chunk_queue, bool real_start=true, u32 offs=0) const;

    //Updates the read chunk
    //only relevent if single_chunk() == true
    //replaces something...
    bool add_next_chunk(ReadBuffer &r);

    #ifdef PYBIND

    #define PY_READ_METH(P) c.def(#P, &ReadBuffer::P);
    #define PY_READ_RPROP(P) c.def_property_readonly(#P, &ReadBuffer::get_##P);
    #define PY_READ_PRM(P) p.def_readwrite(#P, &ReadBuffer::Params::P);

    static void pybind_defs(pybind11::class_<ReadBuffer> &c) {
        PY_READ_METH(empty);
        PY_READ_METH(size); //TODO bind to __len__ ?
        PY_READ_RPROP(id);
        PY_READ_RPROP(start);
        PY_READ_RPROP(end);
        PY_READ_RPROP(duration);
        PY_READ_RPROP(channel);
        PY_READ_RPROP(raw);

        c.def("__len__", &ReadBuffer::size);
        c.def("__getitem__", &ReadBuffer::operator[]);

        pybind11::class_<Params> p(c, "Params");
        PY_READ_PRM(num_channels);
        PY_READ_PRM(bp_per_sec);
        PY_READ_PRM(sample_rate);
        PY_READ_PRM(chunk_time);
        PY_READ_PRM(max_chunks);
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
