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

#ifndef _INCL_FAST5_READER
#define _INCL_FAST5_READER

#include <thread>
#include <vector>
#include <deque>
#include <unordered_set>
#include "read_buffer.hpp"
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#endif

class Fast5Reader {
    public:

    typedef struct {
        std::string fast5_list;
        std::string read_list;
        u32 max_reads, max_buffer;
    } Params;
    static Params const PRMS_DEF;

    typedef struct {
        const char *fast5_list, *read_list, *max_reads, *max_buffer;
    } Docstrs;
    static constexpr Docstrs DOCSTRS = {
        fast5_list : 
            "File containing a list of paths to fast5 files, one per line.",
        read_list : 
            "File containing a list of read IDs. Only these reads will be loaded if specified.",
        max_reads : 
            "Maximum number of reads to load.",
        max_buffer : 
            "Maximum number of reads to store in memory."
    };


    //TODO: remove reduntant constructors

    Fast5Reader();
    Fast5Reader(const Params &p);

    Fast5Reader(u32 max_reads, u32 max_buffer=100);
    Fast5Reader(const std::string &fast5_list, 
                const std::string &read_list="",
                u32 max_reads=0, u32 max_buffer=100);

    void add_fast5(const std::string &fast5_path);

    bool load_fast5_list(const std::string &fname);

    bool add_read(const std::string &read_id);

    bool load_read_list(const std::string &fname);

    ReadBuffer pop_read();
 
    u32 buffer_size();
 
    u32 fill_buffer();
 
    bool all_buffered();
 
    bool empty();

    #ifdef PYBIND

    #define PY_FAST5_METH(N) c.def(#N, &Fast5Reader::N);
    #define PY_FAST5_PRM(P) p.def_readwrite(#P, &Fast5Reader::Params::P);

    static void pybind_defs(pybind11::class_<Fast5Reader> &c) {

        c.def(pybind11::init());
        c.def(pybind11::init<Params>());
        c.def(pybind11::init<u32, u32>());
        c.def(pybind11::init<
                const std::string &, 
                const std::string &, 
                u32, u32>());

        PY_FAST5_METH(add_fast5);
        PY_FAST5_METH(load_fast5_list);
        PY_FAST5_METH(add_read);
        PY_FAST5_METH(load_read_list);
        PY_FAST5_METH(pop_read);
        PY_FAST5_METH(buffer_size);
        PY_FAST5_METH(fill_buffer);
        PY_FAST5_METH(all_buffered);
        PY_FAST5_METH(empty);

        pybind11::class_<Params> p(c, "Params");
        PY_FAST5_PRM(fast5_list);
        PY_FAST5_PRM(read_list);
        PY_FAST5_PRM(max_reads);
        PY_FAST5_PRM(max_buffer);
    }

    #endif

    private:
    Params PRMS;

    enum Format {MULTI, SINGLE, UNKNOWN};
    static const std::string FMT_RAW_PATHS[], FMT_CH_PATHS[];

    bool open_next();

    u32 max_buffer_, total_buffered_, max_reads_;

    std::deque<std::string> fast5_list_;
    std::unordered_set<std::string> read_filter_;

    hdf5_tools::File open_fast5_;
    Format open_fmt_;
    std::deque<std::string> read_paths_;

    std::deque<ReadBuffer> buffered_reads_;
};

#endif
