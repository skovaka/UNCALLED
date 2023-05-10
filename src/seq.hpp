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

#ifndef _INCL_BWAFMI
#define _INCL_BWAFMI

#include <string>
#include <climits>
#include <utility>
#include <cstring>
#include <algorithm>
#include "pore_model.hpp"
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

//From submods/bwa/bwtindex.c
#define BWA_BLOCK_SIZE 10000000

struct SeqRecord {
    std::string name;
    i32 id, length, offset;

    static void pybind(py::module_ &m) {
        py::class_<SeqRecord> c(m, "SeqRecord");
        c.def(py::init<std::string, i32, i32, i32>());
        c.def_readwrite("name", &SeqRecord::name);
        c.def_readwrite("id", &SeqRecord::id);
        c.def_readwrite("length", &SeqRecord::length);
        c.def_readwrite("offset", &SeqRecord::offset);
    }
};

enum class Strand {fwd='+', rev='-', na='*'};
class RefCoord {
    public:
    std::string name;
    std::vector<i64> bounds; 
    //i64 start, end, ref_len;
    Strand strand;
    i32 ref_id;
    //bool fwd;

    RefCoord(std::string _name, std::vector<i64> _bounds, Strand _strand=Strand::na, i32 _ref_id=-1) :
        name(_name),
        bounds(_bounds),
        strand(_strand),
        ref_id(_ref_id) {

        if (bounds.size() == 0 || bounds.size() % 2 != 0) {
            throw std::runtime_error("RefCoord.bounds length must be divisible by 2 and greater than 0");
        }
    }

    RefCoord(std::string _name, IntervalIndex<i64> mpos, bool fwd) :
        name(_name),
        strand(fwd ? Strand::fwd : Strand::rev) {

        if (mpos.get_start() < 0) {
            mpos = mpos.mirror();
        }
        bounds.reserve(mpos.coords.size()*2);

        for (auto b : mpos.coords) {
            bounds.push_back(b.start);
            bounds.push_back(b.end);
        }
    }

    RefCoord(std::string _name, i64 _start, i64 _end, Strand _strand=Strand::na) :
        RefCoord(_name, {_start, _end}, _strand) {}

    RefCoord(std::string _name, std::vector<i64> _bounds, bool fwd) :
        RefCoord(_name, _bounds, fwd ? Strand::fwd : Strand::rev) {}

    RefCoord(std::string _name, i64 _start, i64 _end, bool fwd) :
        RefCoord(_name, {_start,_end}, fwd ? Strand::fwd : Strand::rev) {}

    RefCoord(const RefCoord &rc, bool fwd) :
        RefCoord(rc.name, rc.bounds, fwd ? Strand::fwd : Strand::rev) {}

    RefCoord() : RefCoord("",{-1,-1}) {};

    //RefCoord intersection(RefCoord other) const {
    //    return RefCoord(name, std::max(start, other.start), std::min(end, other.end), strand, ref_id);
    //}

    i64 start() const {
        return bounds.front();
    }

    i64 end() const {
        return bounds.back();
    }

    void set_start(i64 v) {
        bounds[0] = v;
    }

    void set_end(i64 v) {
        bounds[bounds.size()-1] = v;
    }

    bool fwd() const {
        return strand == Strand::fwd;
    }

    bool rev() const {
        return strand == Strand::rev;
    }

    bool stranded() const {
        return strand != Strand::na;
    }

    i64 length() const {
        return end()-start();
    }

    #ifdef PYBIND

    static void pybind_defs(pybind11::module_ m) {
        py::class_<RefCoord> c(m, "_RefCoord");

        c.def(py::init<std::string, IntervalIndex<i64>, bool>());        
        c.def(py::init<std::string, std::vector<i64>>());
        c.def(py::init<std::string, std::vector<i64>, bool>());
        c.def(py::init<std::string, i64, i64>());
        c.def(py::init<std::string, i64, i64, bool>());
        c.def(py::init<const RefCoord &, bool>());
        c.def(py::init<const RefCoord &>());

        c.def("__repr__", 
            [](RefCoord &c) -> std::string {
                std::stringstream ss;
                ss << c.name;
                if (c.start() >= 0) {
                    ss << ":" << c.bounds[0] << "-" << c.bounds[1];
                    for (size_t i = 2; i < c.bounds.size(); i += 2) {
                        ss << ";" << c.bounds[i] << "-" << c.bounds[i+1];
                    }
                }
                return ss.str();
        });


        //c.def("intersection", &RefCoord::intersection);
        c.def("__len__", &RefCoord::length);
        c.def_readwrite("name", &RefCoord::name);
        c.def_readonly("bounds", &RefCoord::bounds);
        c.def_property("start", &RefCoord::start, &RefCoord::set_start);
        c.def_property("end", &RefCoord::end, &RefCoord::set_end);
        //c.def_readwrite("ref_len", &RefCoord::ref_len);
        c.def_property("fwd", &RefCoord::fwd, 
            [](RefCoord &c, bool fwd) -> bool {
                if (fwd) {
                    c.strand = Strand::fwd;
                    return true;
                }
                c.strand = Strand::rev;
                return false;
        });
        c.def_property_readonly("rev", &RefCoord::rev);
        c.def_property_readonly("stranded", &RefCoord::stranded);
    }

    #endif
};


template <typename ModelType>
struct Sequence {//: public DataFrame<typename ModelType::kmer_t, float, u8> {
    using KmerType = typename ModelType::kmer_t;
    //using super = DataFrame<KmerType, u8, float>;

    const ModelType &model;
    const RefCoord coord;
    const KmerType KMER_LEN; //= ModelType::KMER_LEN;
    IntervalIndex<i64> mpos;
    bool is_fwd; //TODO infer from mpos mpos

    //static constexpr typename super::NameArray names = {"ref", "start", "end"}; 
    //typename super::template ColType<0> &kmer = std::get<0>(super::data_);   
    //typename super::template ColType<1> &current = std::get<1>(super::data_);   
    //typename super::template ColType<2> &base = std::get<2>(super::data_);   
    //ValArray<u8> base;

    ValArray<KmerType> kmer; 
    ValArray<float> current;   
    ValArray<bool> splice_mask;   

    Sequence(const ModelType &model_, size_t length) : 
        model(model_), 
        coord("", 0, length),
        KMER_LEN(model_.KMER_LEN),
        mpos({{0,static_cast<i64>(length)}}), 
        is_fwd(true),
        kmer(length), current(length) {}

    //Sequence(const ModelType &model_, IntervalIndex<i64> coords_, bool is_fwd_) : 
    //    model(model_), 
    //    coord("", coords_.get_start(), coords_.get_end(), is_fwd_),
    //    KMER_LEN(model_.KMER_LEN),
    //    mpos(coords_),
    //    is_fwd(is_fwd_), 
    //    kmer(mpos.length), 
    //    current(mpos.length) {
    //    init_current();
    //}

    Sequence(const ModelType &model_, const std::string &seq, RefCoord &coords_) : 
            model(model_), 
            coord(coords_),
            KMER_LEN(model_.KMER_LEN),
            is_fwd(coords_.fwd()) {

        auto rev = is_fwd == model.PRMS.reverse;
        auto comp = !is_fwd;

        if (!rev) {
            auto start = coords_.bounds[0] + model_.PRMS.shift;
            for (size_t i = 1; i < coords_.bounds.size()-1; i += 2) {
                auto end = coords_.bounds[i];
                mpos.append(start, end);
                start = coords_.bounds[i+1];
            }
            auto end = coords_.bounds.back() - KMER_LEN + model_.PRMS.shift + 1;
            mpos.append(start, end);
        } else {
            auto start = -coords_.bounds.back() + model_.PRMS.shift;
            for (size_t i = coords_.bounds.size()-2; i > 0; i -= 2) {
                auto end = -coords_.bounds[i];
                mpos.append(start, end);
                start = -coords_.bounds[i-1];
            }
            auto end = -coords_.bounds[0] - KMER_LEN + model_.PRMS.shift + 1;
            mpos.append(start, end);
        }

        if (mpos.coords.size() > 1) {
            splice_mask = ValArray<bool>(true, mpos.length);
            size_t i = 0;
            for (size_t j = 0; j < mpos.coords.size()-1; j++) {
                i += mpos.coords[j].length();
                for (int sh = -model_.PRMS.shift; sh < model_.PRMS.shift; sh++) {
                    splice_mask[i+sh] = false;
                }
            }
        }

        if (mpos.length != seq.size() - KMER_LEN + 1) {
            throw std::runtime_error("Size mismatch");
        }
        current = ValArray<float>(mpos.length);
        kmer = model.str_to_kmers(seq, rev, comp); 
        init_current();
    }

    Sequence(const ModelType &model_, const std::string &seq) :
            Sequence(model_, seq.size()-KMER_LEN+1) {
        kmer = model.str_to_kmers(seq);
        init_current();
    }

    //Sequence(const ModelType &model_, u8 *seq, size_t start, size_t end) :
    //        Sequence(model_, end-start-KMER_LEN+1) {
    //    kmer = model.pacseq_to_kmers(seq, start, end);
    //    init_current();
    //}
    //
    bool is_spliced() const {
        return splice_mask.size() > 0;
    }

    KmerType get_kmer(i64 r) const {
        auto i = mpos.get_index(r);
        return kmer[i];
    }

    KmerType get_current(i64 r) const {
        auto i = mpos.get_index(r);
        return current[i];
    }

    void init_current() {
        for (size_t i = 0; i < size(); i++) {
            current[i] = model.kmer_current(kmer[i]);
        }
    }

    size_t size() const {
        return kmer.size();
    }

    const ModelType &get_model() {
        return model;
    }

    static void pybind(py::module &m, std::string suffix) {
        py::class_<Sequence> c(m, ("Sequence"+suffix).c_str());
        //auto c = super::template pybind<Sequence>(m, ("Sequence"+suffix).c_str(), false);

        c.def(py::init<const ModelType &, const std::string &>());
        c.def(py::init<const ModelType &, const std::string &, RefCoord &>());
        c.def("__len__", &Sequence::size);
        c.def_property_readonly("model", &Sequence::get_model);
        c.def_readonly("K", &Sequence::KMER_LEN);
        c.def_readonly("mpos", &Sequence::mpos);
        //c.def_readonly("ref_id", &Sequence::ref_id);
        c.def_readonly("kmer", &Sequence::kmer);
        c.def_readonly("current", &Sequence::current);
        c.def_readonly("is_fwd", &Sequence::is_fwd);
        c.def_readonly("coord", &Sequence::coord);
        c.def_readonly("splice_mask", &Sequence::splice_mask);
        //c.def("kmer_to_str", py::vectorize(&Sequence::kmer_to_str));
        c.def("is_spliced", &Sequence::is_spliced);
        c.def("get_kmer", py::vectorize(&Sequence::get_kmer));
        c.def("get_current", py::vectorize(&Sequence::get_current));
    }
};


#endif
