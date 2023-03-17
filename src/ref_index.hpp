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
#include <bwa/bwa.h>
#include <bwa/utils.h>
#include <pdqsort.h>
#include <zlib.h>
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

enum class Strand {fwd='+', rev='-', na='*'};

class RefCoord {
    public:
    std::string name;
    i64 start, end, ref_len;
    Strand strand;
    i32 ref_id;
    //bool fwd;
    
    RefCoord(std::string _name, i64 _start, i64 _end, Strand _strand=Strand::na, 
              i32 _ref_id=-1, i64 _ref_len=-1) :
        name(_name),
        start(_start),
        end(_end),
        ref_len(_ref_len),
        strand(_strand),
        ref_id(_ref_id) {}

    RefCoord(std::string _name, i64 _start, i64 _end, bool fwd) :
        RefCoord(_name,_start,_end, fwd ? Strand::fwd : Strand::rev) {}

    RefCoord(const RefCoord &rc, bool fwd) :
        RefCoord(rc.name,rc.start,rc.end, fwd ? Strand::fwd : Strand::rev) {}

    RefCoord() : RefCoord("",-1,-1) {};

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
        return end-start;
    }

    #ifdef PYBIND

    static void pybind_defs(pybind11::module_ m) {
        py::class_<RefCoord> c(m, "_RefCoord");

        c.def(py::init<std::string, i64, i64>());
        c.def(py::init<std::string, i64, i64, bool>());
        c.def(py::init<const RefCoord &, bool>());
        c.def(py::init<const RefCoord &>());

        c.def("__repr__", 
            [](RefCoord &c) -> std::string {
                if (c.start < 0) {
                    return c.name;
                }
                auto s = c.name + ":" + 
                         std::to_string(c.start) + "-" +
                         std::to_string(c.end);
                if (c.stranded()) {
                    return s + ":" + std::string(1, static_cast<char>(c.strand));
                }
                return s;
        });


        c.def_readwrite("name", &RefCoord::name);
        c.def_readwrite("start", &RefCoord::start);
        c.def_readwrite("end", &RefCoord::end);
        c.def_readwrite("ref_len", &RefCoord::ref_len);
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

template<class ModelType>
class RefIndex {
    public:

    using KmerType = typename ModelType::kmer_t;
    ModelType model_;
    const KmerLen K;

    using Range = std::pair<u64, u64>;

    RefIndex(ModelType &model) :
        model_(model),
        K(model.KMER_LEN),
        index_(NULL),
        bns_(NULL),
        pacseq_(NULL),
        kmer_ranges_(model_.KMER_COUNT),
        loaded_(false),
        size_(0) {}

    RefIndex(ModelType &model, const std::string &prefix, bool pacseq=false, bool bwt=false) : 
            RefIndex(model) {
        if (!prefix.empty()) {
            bns_ = bns_restore(prefix.c_str());
            size_ = 2 * (bns_->l_pac);
            if (bwt) load_index(prefix);
            if (pacseq) load_pacseq();
        }
    }

    void load_index(const std::string &prefix) {
        std::string bwt_fname = prefix + ".bwt",
                    sa_fname = prefix + ".sa";

        if (bns_ == NULL) {
            bns_ = bns_restore(prefix.c_str());
        }

        index_ = bwt_restore_bwt(bwt_fname.c_str());
        bwt_restore_sa(sa_fname.c_str(), index_);

        for (KmerType k = 0; k < kmer_ranges_.size(); k++) {

            auto r = get_base_range(model_.kmer_head(k));
            for (u8 i = 1; i < K; i++) {
                r = get_neighbor(r, model_.kmer_base(k, i));
            }

            kmer_ranges_[k] = r;
        }

        loaded_ = true;
    }

    static void create(const std::string &fasta_fname, 
                       const std::string &prefix = "",
                       bool no_bwt=false) {

        std::string prefix_auto = prefix.empty() ? fasta_fname : prefix;

        if (no_bwt) {
            gzFile fp = xzopen(fasta_fname.c_str(), "r");
            bns_fasta2bntseq(fp, prefix_auto.c_str(), 0);
            err_gzclose(fp);

        } else {
            bwa_idx_build(fasta_fname.c_str(), 
                          prefix.c_str(), 
                          BWTALGO_AUTO,
                          BWA_BLOCK_SIZE);
        }
    }

    bool bwt_loaded() {
        return loaded_;
    }

    void load_pacseq() {
        if (!pacseq_loaded()) {
            //Copied from bwa/bwase.c
            pacseq_ = (u8*) calloc(bns_->l_pac/4+1, 1);
            err_fread_noeof(pacseq_, 1, bns_->l_pac/4+1, bns_->fp_pac);
        }   
    }

    void destroy() {
        if (index_ != NULL) { 
            bwt_destroy(index_);
        }
        if (bns_ != NULL) { 
            bns_destroy(bns_);
        }
        if (pacseq_loaded()) {
            free(pacseq_);
        }
    }

    Range get_neighbor(Range r1, u8 base) const {
        u64 os, oe;
        bwt_2occ(index_, r1.first - 1, r1.second, base, &os, &oe);
        return Range(index_->L2[base] + os + 1, index_->L2[base] + oe);
    }

    Range get_kmer_range(KmerType kmer) const {
        return kmer_ranges_[kmer];
    }

    i64 get_kmer_count(KmerType kmer) const {
        auto r = kmer_ranges_[kmer];
        return r.first - r.second + 1;
    }

    Range get_base_range(u8 base) const {
        return {index_->L2[base], index_->L2[base+1]};
    }

    i64 fm_to_mpos(i64 fm) {
        //return size() - fm_to_pac(fm);
        return -fm_to_pac(fm)-1;
    }

    i64 fm_to_pac(i64 fm) const {
        return bwt_sa(index_, fm);
    }

    size_t size() const {
        //return index_->seq_len;
        return size_;
    }

    i32 get_ref_id(const std::string &ref_name) {
        for (int i = 0; i < bns_->n_seqs; i++) {
            if (strcmp(bns_->anns[i].name, ref_name.c_str()) == 0) {
                return i;
            }
        }
        return -1;
    }


    i32 pac_to_ref_id(i64 pac) {
        return bns_pos2rid(bns_, pac);
    }

    bool is_mpos_flipped(i64 i) const {
        return i < 0;
    }

    bool is_mpos_fwd(i64 i, bool is_rna) const {
        return is_mpos_flipped(i) == is_rna;
    }

    i64 mpos_to_pac(i32 ref_id, i64 mpos) {
        if (is_mpos_flipped(mpos)) {
            mpos = -mpos-1;
        }
        return get_pac_offset(ref_id) + mpos;
    }

    std::string get_ref_name(u32 rid) {
        return bns_->anns[rid].name;
    }


    i64 get_ref_len(u32 rid) const {
        return bns_->anns[rid].len;
    }

    i64 get_sa_loc(const std::string &name, i64 coord) {
        for (int i = 0; i < bns_->n_seqs; i++) {
            if (strcmp(bns_->anns[i].name, name.c_str()) == 0) {
                return bns_->anns[i].offset + coord;
            }
        }
        return 0;
    }


    RefCoord mpos_to_pos_coord(i64 mpos_start, i64 mpos_end, bool read_fwd) {
        bool flip = mpos_start < 0;// static_cast<i32>(size() / 2);

        i64 pac_st; //TODO rename
        if (is_mpos_flipped(mpos_start)) { 
            pac_st = mpos_to_pac(0, mpos_end)+1;
        } else {
            pac_st = mpos_start;
        }

        i32 rid = bns_pos2rid(bns_, pac_st);
        if (rid < 0) return RefCoord();
        //assert(rid >= 0);

        auto name = std::string(bns_->anns[rid].name);
        auto start = pac_st - bns_->anns[rid].offset;
        auto end = start + (mpos_end-mpos_start);
        auto strand = flip != read_fwd ? Strand::fwd : Strand::rev;

        return RefCoord(name,start,end,strand,rid,bns_->anns[rid].len);
    }

    std::vector< std::pair<std::string, i64> > get_seqs() const {
        std::vector< std::pair<std::string, i64> > seqs;

        for (i32 i = 0; i < bns_->n_seqs; i++) {
            bntann1_t ann = bns_->anns[i];
            std::string name = std::string(ann.name);
            seqs.push_back( std::pair<std::string, i64>(name, ann.len) );
        }

        return seqs;
    }

    i64 pos_to_pac(std::string name, i64 coord) {
        i32 i;
        for (i = 0; i < bns_->n_seqs; i++) {
            if (strcmp(name.c_str(), bns_->anns[i].name) == 0)
                return bns_->anns[i].offset + coord;
        }
        return INT_MAX;
    } 

    //TODO overload same function below
    i64 get_pac_offset(const std::string &ref_name) {
        return get_pac_offset(get_ref_id(ref_name));
    }

    i64 get_pac_offset(i32 ref_id) {
        return bns_->anns[ref_id].offset;
    }

    bool pacseq_loaded() const {
        return pacseq_ != NULL;
    }

    std::pair<i64, i64> pos_to_mpos(i64 st, i64 en, bool is_fwd, bool is_rna) {
        auto flip = is_fwd == is_rna;
        if (!flip) return {st, en};
        return {-en, -st};
    }

    i64 pos_to_mpos(i64 ref, bool is_fwd, bool is_rna) {
        if (is_fwd == is_rna) return -ref-1;
        return ref;
    }

    i64 mpos_to_pos(i64 mpos) {
        if (mpos < 0) return -mpos-1;
        return mpos;
    }

    i64 pac_to_pos(i64 pac) {
        i32 rid = bns_pos2rid(bns_, pac);
        return pac - bns_->anns[rid].offset;
    }

    u8 get_base(i64 pac, bool comp=false) {
        if (pac < 0 || pac > static_cast<i64>(size() / 2)) { //TODO better size definition
            throw std::out_of_range("Base out of range: " + std::to_string(pac));
        }
        size_t i = pac >> 2,
               shift = ((pac & 3) ^ 3) << 1;
        return ((pacseq_[i] >> shift) & 3) ^ (comp | (comp<<1));
    }
    
    private:
    KmerType next_kmer(KmerType kmer, i64 pac, bool comp) {
        auto b = get_base(pac, comp);
        auto ret = model_.kmer_neighbor(kmer, b);
        return ret;
    }

    void next_kmer(std::vector<KmerType> &kmers, i64 pac, bool comp) {
        kmers.push_back(next_kmer(kmers.back(), pac, comp));
    }

    public:
    KmerType get_kmer(i64 pac, bool comp) {
        KmerType kmer = 0;
        for (auto i = pac; i < pac + K; i++) {
            kmer = next_kmer(kmer, i, comp);
        }
        return kmer;
    }

    //std::vector<KmerType> get_kmers(std::vector<std::pair<i64, i64>> pac_blocks, bool rev, bool comp) {
    //    size_t len = 0;
    //    for (auto &b : pac_blocks) {
    //        len += b.second - b.first;
    //    }
    //    len -= K - 1;
    //    std::vector<KmerType> kmers(len);
    //    for (auto &b : pac_blocks) {
    //        get_kmers(b.first, b.second, rev, comp, kmers);
    //    }
    //    return kmers;
    //}

    std::vector<KmerType> get_kmers(i64 pac_start, i64 pac_end, bool rev, bool comp) {
        std::vector<KmerType> kmers(pac_end-pac_start-K+1);
        get_kmers(pac_start, pac_end, rev, comp, kmers);
        return kmers;
    }

    template <typename Container>
    size_t get_kmers(i64 pac_start, i64 pac_end, bool rev, bool comp, Container &kmers, size_t k=0) {
        if (!rev) {
            auto i = pac_start;
            if (k == 0) {
                kmers[k++] = get_kmer(i, comp);
                //kmers.push_back(get_kmer(i, comp));
                i += K;
            }
            for (; i < pac_end; i++) {
                kmers[k++] = next_kmer(kmers[k-1], i, comp);
                //next_kmer(kmers, i, comp);
            }
        } else {
            auto i = pac_end-1;
            if (k == 0) {
                kmers[k++] = model_.kmer_rev(get_kmer(pac_end-K, comp));
                i -= K;
            }
            //for (auto i = pac_end-K-1; i >= pac_start; i--) {
            for (; i >= pac_start; i--) {
                kmers[k++] = next_kmer(kmers[k-1], i, comp);
                //next_kmer(kmers, i, comp);
            }
        }
        return k;
    }

    std::vector<KmerType> get_kmers(const std::string &name, i64 start, i64 end, bool rev=false, bool comp=false) {
        i64 pac_start = pos_to_pac(name, start),
            pac_end = pos_to_pac(name, end);
        return get_kmers(pac_start, pac_end, rev, comp);
    }

    //std::vector<KmerType> get_kmers(i32 ref_id, std::vector<std::pair<i64, i64>> mpos_blocks, bool is_rna) {
    //    size_t len = 0;
    //    for (auto &b : mpos_blocks) {
    //        len += b.second - b.first;
    //    }
    //    len -= K - 1;
    //    std::vector<KmerType> kmers(len);
    //    size_t i = 0;
    //    for (auto &b : mpos_blocks) {
    //        i = get_kmers(ref_id, b.first, b.second, is_rna, kmers, i);
    //    }
    //    return kmers;
    //}

    Sequence<ModelType> get_kmers(ModelType &model, i32 ref_id, IntervalIndex<i64> &mpos_blocks, bool is_rna) {
        auto mpos_end = mpos_blocks.coords.back().end-1;  
        bool rev = is_mpos_flipped(mpos_end);
        auto seq_idx = mpos_blocks;//.islice(trim_st, mpos_blocks.length-trim_en);
        
        Sequence<ModelType> seq(model, ref_id, seq_idx, rev == is_rna);
        auto lpad = (K - 1) / 2, rpad = K - lpad - 1;
        size_t k = 0;
        for (size_t i = 0; i < mpos_blocks.interval_count()-1; i++) {
            auto &c = mpos_blocks.coords[i];
            k = get_kmers(ref_id, c.start-lpad, c.end, is_rna, seq.kmer, k);
            lpad = 0;
        }
        auto &c = mpos_blocks.coords.back();
        get_kmers(ref_id, c.start-lpad, c.end+rpad, is_rna, seq.kmer, k);

        seq.init_current();
        return seq;
    }
    
    std::vector<KmerType> get_kmers(i32 ref_id, i64 mpos_start, i64 mpos_end, bool is_rna) {
        std::vector<KmerType> kmers(mpos_end-mpos_start - K + 1);
        get_kmers(ref_id, mpos_start, mpos_end, is_rna, kmers, 0);
        return kmers;
    }

    template <typename Container>
    size_t get_kmers(i32 ref_id, i64 mpos_start, i64 mpos_end, bool is_rna, Container &kmers, size_t k) {
        bool rev = is_mpos_flipped(mpos_end-1);
        bool comp = rev != is_rna;
        i64 pac_start = mpos_to_pac(ref_id, mpos_start),
            pac_end = mpos_to_pac(ref_id, mpos_end);

        if (rev) {
            return get_kmers(pac_end+1, pac_start+1, true, comp, kmers, k);
        } else {
            return get_kmers(pac_start, pac_end, false, comp, kmers, k);
        }
    }

    using FwdRevCoords = std::pair< std::vector<i64>, std::vector<i64> >;


    #ifdef PYBIND

    #define PY_BWA_INDEX_METH(P) c.def(#P, &RefIndex<ModelType>::P);
    #define PY_BWA_INDEX_VEC(P) c.def(#P, py::vectorize(&RefIndex<ModelType>::P));

    static void pybind_defs(pybind11::module_ m, const std::string &suffix) {
        py::class_<RefIndex<ModelType>> c(m, ("RefIndex" + suffix).c_str());

        c.def(py::init<ModelType &>());
        c.def(py::init<ModelType &, const std::string &>());
        c.def(py::init<ModelType &, const std::string &, bool>());
        c.def(py::init<ModelType &, const std::string &, bool, bool>());
        PY_BWA_INDEX_METH(create);
        PY_BWA_INDEX_METH(load_index);
        PY_BWA_INDEX_METH(bwt_loaded);
        PY_BWA_INDEX_METH(load_pacseq);
        PY_BWA_INDEX_METH(destroy);
        PY_BWA_INDEX_METH(get_neighbor);
        PY_BWA_INDEX_METH(get_kmer_range);
        PY_BWA_INDEX_METH(get_kmer_count);
        c.def("get_kmer_count", py::vectorize(&RefIndex<ModelType>::get_kmer_count));
        PY_BWA_INDEX_METH(get_base_range);
        PY_BWA_INDEX_METH(fm_to_pac);
        PY_BWA_INDEX_METH(fm_to_mpos);
        PY_BWA_INDEX_VEC(mpos_to_pac);
        PY_BWA_INDEX_VEC(pac_to_pos);
        PY_BWA_INDEX_METH(pos_to_pac);
        PY_BWA_INDEX_METH(size);
        PY_BWA_INDEX_METH(mpos_to_pos_coord);
        PY_BWA_INDEX_METH(get_seqs);
        PY_BWA_INDEX_METH(pacseq_loaded);
        PY_BWA_INDEX_METH(get_sa_loc);
        PY_BWA_INDEX_METH(get_ref_name);
        PY_BWA_INDEX_METH(get_ref_len);
        c.def("get_pac_offset", static_cast<i64 (RefIndex::*)(i32)> (&RefIndex::get_pac_offset) );
        c.def("get_pac_offset", static_cast<i64 (RefIndex::*)(const std::string &)> (&RefIndex::get_pac_offset) );
        PY_BWA_INDEX_METH(is_mpos_fwd);
        PY_BWA_INDEX_METH(is_mpos_flipped);
        PY_BWA_INDEX_VEC(get_base);
        c.def("pac_to_ref_id", static_cast<i32 (RefIndex::*)(i64)> (&RefIndex::pac_to_ref_id) );
        c.def("get_ref_id", static_cast<i32 (RefIndex::*)(const std::string &)> (&RefIndex::get_ref_id));
        //c.def("get_kmers_new", &RefIndex::get_kmers_new);
        c.def("pos_to_mpos", static_cast<std::pair<i64, i64> (RefIndex::*)(i64, i64, bool, bool)> (&RefIndex::pos_to_mpos));
        c.def("pos_to_mpos", pybind11::vectorize(static_cast<i64 (RefIndex::*)(i64, bool, bool)> (&RefIndex::pos_to_mpos)));
        c.def("mpos_to_pos", pybind11::vectorize(&RefIndex::mpos_to_pos));
        //c.def("get_kmers", static_cast< std::vector<KmerType> (RefIndex::*)(i64, i64)> (&RefIndex::get_kmers) );
        
        c.def("get_kmers", 
            static_cast< std::vector<KmerType> (RefIndex::*)(i64, i64, bool, bool)> (&RefIndex::get_kmers),
            py::arg("pac_start"), py::arg("pac_end"), py::arg("rev"), py::arg("comp"));

        //c.def("get_kmers", 
        //    static_cast< std::vector<KmerType> (RefIndex::*)(std::vector<std::pair<i64, i64>>, bool, bool)> (&RefIndex::get_kmers),
        //    py::arg("pac_blocks"), py::arg("rev"), py::arg("comp"));

        //c.def("get_kmers", 
        //    static_cast< std::vector<KmerType> (RefIndex::*)(i32, std::vector<std::pair<i64, i64>>, bool)> (&RefIndex::get_kmers),
        //    py::arg("ref_id"), py::arg("mpos_blocks"), py::arg("is_rna"));

        c.def("get_kmers", 
            static_cast< Sequence<ModelType> (RefIndex::*)(ModelType&, i32, IntervalIndex<i64>&, bool)> (&RefIndex::get_kmers),
            py::arg("model"), py::arg("ref_id"), py::arg("mpos_blocks"), py::arg("is_rna"));

        c.def("get_kmers", 
            static_cast< std::vector<KmerType> (RefIndex::*)(const std::string &, i64, i64, bool, bool)> (&RefIndex::get_kmers),

            py::arg("name"), py::arg("start"), py::arg("end"), py::arg("rev")=false, py::arg("comp")=false);
        c.def("get_kmers", 
            static_cast< std::vector<KmerType> (RefIndex::*)(i32, i64, i64, bool)> (&RefIndex::get_kmers),
            py::arg("ref_id"), py::arg("miref_start"), py::arg("miref_end"), py::arg("is_rna"));

    }

    #endif

    private:
    bwt_t *index_;
    bntseq_t *bns_;
    u8 *pacseq_;
    std::vector<Range> kmer_ranges_;
    bool loaded_;
    size_t size_;
};


#endif
