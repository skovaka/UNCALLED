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

#include <string>
#include <vector>
#include <unordered_set>
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

class ReadBuffer {
    public:

    //typedef struct {
    //    u16 num_channels;
    //    float bp_per_sec;
    //    float sample_rate;
    //    bool seq_fwd; //If true reads start at 5', end at 3' (DNA)
    //} Params;

    //static Params PRMS;

    ReadBuffer();

    ReadBuffer(const std::string &id, u16 channel, u32 number, u64 start_time, 
          const std::string &dtype, const std::string &raw_str);
    ReadBuffer(const std::string &id, u16 channel, u32 number, u64 start_time, 
          const std::vector<float> &raw_data, u32 raw_st, u32 raw_len);
    
    #ifdef PYBIND
    ReadBuffer(const std::string &id, u16 channel, u32 number, u64 start_time, const py::array_t<float> &raw_data);
    #endif

    //Returns start sample of read
    u64 get_start() const;

    //Returns last known sample of read
    u64 get_end() const;

    //Returns size of buffered signal in samples
    u64 size() const;

    void set_signal(const float *ptr, size_t len);
    
    //Returns sample at specified index
    float &operator[](u32 i);

    //Clears the signal buffer
    void clear();
    
    //Returns true if no signal is in buffer
    bool empty() const;

    template <typename Container>
    void set_moves(const Container &moves_, u32 template_start_, u32 stride) {
        moves = ValArray<bool>(moves_.begin(), moves_.size());
        template_start = template_start_;
        move_stride = stride;
        bc_loaded = true;
    }

    #ifdef PYBIND

    #define PY_READ_METH(P, D) c.def(#P, &ReadBuffer::P, D);
    #define PY_READ_RPROP(P, D) c.def_property_readonly(#P, &ReadBuffer::get_##P, D);
    #define PY_READ_PRM(P, D) p.def_readwrite(#P, &ReadBuffer::Params::P, D);

    #define PY_READ_ATTR(P,D) c.def_readonly(#P, &ReadBuffer::P, D);

    void set_moves(py::array_t<bool> moves_py, u32 template_start, u32 stride) {
        PyArray<bool> moves_(moves_py);
        set_moves(moves_, template_start, stride);
    }

    //static void pybind_defs(pybind11::class_<ReadBuffer> &c) {
    static void pybind_defs(pybind11::module_ m) {
        py::class_<ReadBuffer> c(m, "ReadBuffer");
        c.def(pybind11::init<ReadBuffer>());
        c.def(pybind11::init<const std::string &, u16, u32, u64, const py::array_t<float> &>());
    
        c.def_readwrite("filename", &ReadBuffer::filename);
        PY_READ_METH(empty, "Returns true if read has no data");
        PY_READ_RPROP(end, "Read end sample relative to start of the run");
        PY_READ_RPROP(start, "Read start sample relative to start of the run");
        PY_READ_ATTR(id, "Read ID (name)");
        PY_READ_ATTR(channel, "Channel where the read was sequenced");
        PY_READ_ATTR(number, "Read number");
        PY_READ_ATTR(signal, "Read signal");
        //c.def_property_readonly("signal", [](ReadBuffer &r) -> pybind11::array_t<float> {
        //     return pybind11::array_t<float>(r.signal.size(), r.signal.data());
        //}, "Read Signal");


        c.def("set_moves", 
            static_cast<void (ReadBuffer::*) (py::array_t<bool>,u32,u32)>(&ReadBuffer::set_moves), "Set move data");
        PY_READ_ATTR(bc_loaded, "True if basecalling data loaded");
        PY_READ_ATTR(moves, "Guppy BC event moves");
        PY_READ_ATTR(move_stride, "Guppy BC event length");
        PY_READ_ATTR(template_start, "Sample where guppy basecalling starts");

        c.def("__len__", &ReadBuffer::size);
        c.def("__getitem__", &ReadBuffer::operator[]);

        //pybind11::class_<Params> p(c, "Params");
        //PY_READ_PRM(flowcell, "Flowcell ID");
        //PY_READ_PRM(kit, "Kit ID");
        //PY_READ_PRM(num_channels, "Number of channels on device");
        //PY_READ_PRM(bp_per_sec, "Expected bases sequenced per second");
        //PY_READ_PRM(sample_rate, "Number of raw samples per second");
        //PY_READ_PRM(seq_fwd, "If true indicates that the signal moves 5' -> 3'");
    }

    #endif

    //Source source_;
    std::string id, filename;
    u16 channel;
    u32 number;
    u64 start_sample; //TODO no raw len
    ValArray<float> signal; //TODO store one

    //IntervalIndex<i32> events;
    //ValArray<float> event_mean, event_stdv;

    bool bc_loaded = false;
    u32 template_start, move_stride;
    ValArray<bool> moves;

    friend bool operator< (const ReadBuffer &r1, const ReadBuffer &r2);
};

bool operator< (const ReadBuffer &r1, const ReadBuffer &r2);

#endif
