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
#include <unordered_map>
#include "fast5_read.hpp"
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

/* Fast5Iter and Fast5Dict extend Fast5Reader
 * 

*/

class Fast5Reader {
    public:

    struct Params {
        std::vector<std::string> fast5_files, read_filter;
        std::string fast5_index;
        u32 max_reads, max_buffer;
        bool recursive, load_bc;

        bool load_fast5_files(const std::string &filename) {
            fast5_files = read_txt_file(filename);
            return !fast5_files.empty();
        }

        bool load_read_filter(const std::string &filename) {
            read_filter = read_txt_file(filename);
            return !read_filter.empty();
        }
    };
    static Params const PRMS_DEF;

    Fast5Reader();
    Fast5Reader(const Params &p);

    u32 add_fast5(const std::string &fast5_path);

    bool load_fast5_files(const std::string &fname);


    protected:
    Params PRMS;

    enum class Format {MULTI=0, SINGLE=1, UNKNOWN};

    bool open_fast5(u32 i);

    std::string get_single_raw_path();
    std::string get_read_id(const std::string &raw_path);
    Fast5Read::Paths get_subpaths(const std::string &base);

    static std::vector<std::string> read_txt_file(const std::string &filename) {
        std::vector<std::string> ret;

        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: failed to open \"" << filename << "\".\n";
            return ret;
        }

        std::string line;
        while (getline(file, line)) {
            ret.push_back(line);
        }

        return ret;
    }

    std::vector<std::string> fast5_files_;
    hdf5_tools::File fast5_file_; //rename fast5_file_
    std::vector<std::string> root_ls_;
    Format fmt_; //rename file_fmt_
    u32 fast5_idx_;

    std::deque<std::string> read_paths_;

    std::deque<Fast5Read> buffered_reads_;

    #ifdef PYBIND

    #define PY_FAST5_READER_METH(N,D) c.def(#N, &Fast5Reader::N, D);
    #define PY_FAST5_READER_PRM(P, D) p.def_readwrite(#P, &Fast5Reader::Params::P, D);

    public:

    static void pybind_defs(pybind11::class_<Fast5Reader> &c) {

        PY_FAST5_READER_METH(add_fast5, "Adds fast5 filename");

        pybind11::class_<Params> p(c, "Params");
        //p.def_readwrite(PARAM_META[0], &Fast5Reader::Params::P);

        //DPRM(p, "fast5_files", &Fast5Reader::Params::fast5_files);
        PY_FAST5_READER_PRM(fast5_files, 
                "List of paths to any combination of:\n"
                " 1. fast5 files\n"
                " 2. directories to search for fast5 files (optionally recursive)\n"
                " 3. text file(s) listing one fast5 file or directory to search per line"
        );
                
        PY_FAST5_READER_PRM(read_filter, "List of read IDs and/or text file(s) containing one read ID per line");
        PY_FAST5_READER_PRM(fast5_index, 
            "Filename mapping between read IDs and fast5 files \n"
            "Can be sequencing summary output by basecaller, "
            "\"filename_mapping.txt\" from ont_fast5_api, or "
            "nanopolish \"*.index.readdb\" file"
        );
        PY_FAST5_READER_PRM(max_reads, "Maximum number of reads to load.");
        PY_FAST5_READER_PRM(max_buffer, "Maximum number of reads to store in memory.");
        PY_FAST5_READER_PRM(load_bc, "Will load basecaller data if true");
        PY_FAST5_READER_PRM(recursive, "Will search directories recursively for any file ending in \".fast5\" if true");
            
    }

    #endif
};

class Fast5Iter : public Fast5Reader {
    public:
    Fast5Iter();
    Fast5Iter(const Params &p);
    Fast5Iter(const std::vector<std::string> &fast5s, const std::vector<std::string> &reads = {}, const Params &p=PRMS_DEF);
    Fast5Iter(const std::string &fast5_list, 
                const std::string &read_filter="");

    bool add_read(const std::string &read_id);

    bool load_read_filter(const std::string &fname);

    Fast5Read next_read();
 
    bool empty();

    private:

    bool open_fast5(u32 i);
 
    bool all_buffered();
 
    u32 fill_buffer();
 
    u32 buffer_size();

    u32 max_buffer_, 
        total_buffered_, 
        max_reads_;

    std::unordered_set<std::string> read_filter_;

    #ifdef PYBIND

    #define PY_FAST5_ITER_METH(N,D) c.def(#N, &Fast5Iter::N, D);
    #define PY_FAST5_ITER_PRM(P, D) p.def_readwrite(#P, &Fast5Iter::Params::P, D);

    public:

    static void pybind_defs(pybind11::class_<Fast5Iter, Fast5Reader> &c) {

        c.def(pybind11::init());
        c.def(pybind11::init<Params>());
        c.def(pybind11::init<
                const std::vector<std::string> &,
                const std::vector<std::string> &,
                const Params &>());
        c.def(pybind11::init<
                const std::vector<std::string> &, 
                const std::vector<std::string> &
        >());
        c.def(pybind11::init<
                const std::vector<std::string> &
        >());

        
        PY_FAST5_ITER_METH(add_read, "Adds a read ID to the read filter");
        PY_FAST5_ITER_METH(load_read_filter, "Loads a list of read IDs from a text file to add to the read filter");
        PY_FAST5_ITER_METH(next_read, "Returns the next read");
        PY_FAST5_ITER_METH(empty, "Returns True if all fast5 reads have been output");
        PY_FAST5_ITER_METH(buffer_size, "Returns the number of fast5 reads stored in the buffer");
    }
    #endif
    
};

class Fast5Dict : public Fast5Reader {
    public:

    using Fast5ReadMap = 
        std::unordered_map<std::string, std::vector<std::string>>;

    Fast5Dict();
    Fast5Dict(Fast5ReadMap fast5_map, const Params &p=PRMS_DEF);
    Fast5Dict(const Params &p);

    void add_read(const std::string &read_id, u32 fast5_idx); 
    bool load_index(const std::string &filename);

    Fast5Read operator[](const std::string &read_id);

    std::string get_read_file(const std::string &read_id);

    private:
    std::unordered_map<std::string, u32> reads_;

    #ifdef PYBIND

    #define PY_FAST5_DICT_METH(N,D) c.def(#N, &Fast5Dict::N, D);
    #define PY_FAST5_DICT_PRM(P, D) p.def_readwrite(#P, &Fast5Dict::Params::P, D);

    public:

    static void pybind_defs(pybind11::class_<Fast5Dict, Fast5Reader> &c) {

        c.def(pybind11::init());
        c.def(pybind11::init<Fast5ReadMap, Params>());
        c.def(pybind11::init<Fast5ReadMap>());
        
        PY_FAST5_DICT_METH(load_index, "");
        PY_FAST5_DICT_METH(add_read, "");
        PY_FAST5_DICT_METH(get_read_file, "");
        //c.def("add_read", pybind11::vectorize(&Fast5Dict::add_read), "");
        c.def("__getitem__", &Fast5Dict::operator[]);
    }
    #endif
};

#endif
