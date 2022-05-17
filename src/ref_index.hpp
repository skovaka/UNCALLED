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
#include "nt.hpp"
#include "range.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
#endif

//From submods/bwa/bwtindex.c
#define BWA_BLOCK_SIZE 10000000

struct RefLoc {
    i32 ref_id;
    std::string ref_name;
    i64 ref_len, start, end;
    bool fwd;
};

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
    static constexpr auto K = ModelType::KMER_LEN;

    RefIndex() :
        index_(NULL),
        bns_(NULL),
        pacseq_(NULL),
        kmer_ranges_(ModelType::KMER_COUNT),
        loaded_(false),
        size_(0) {}

    RefIndex(const std::string &prefix, bool pacseq=false, bool bwt=true) : 
            index_(NULL),
            bns_(NULL),
            pacseq_(NULL),
            kmer_ranges_(ModelType::KMER_COUNT),
            loaded_(false),
            size_(0) {
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

            Range r = get_base_range(ModelType::kmer_head(k));
            for (u8 i = 1; i < K; i++) {
                r = get_neighbor(r, ModelType::kmer_base(k, i));
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
        bwt_2occ(index_, r1.start_ - 1, r1.end_, base, &os, &oe);
        return Range(index_->L2[base] + os + 1, index_->L2[base] + oe);
    }

    Range get_kmer_range(KmerType kmer) const {
        return kmer_ranges_[kmer];
    }

    i64 get_kmer_count(KmerType kmer) const {
        return kmer_ranges_[kmer].length();
    }

    Range get_base_range(u8 base) const {
        return Range(index_->L2[base], index_->L2[base+1]);
    }

    i64 fm_to_mref(i64 fm) {
        return size() - fm_to_pac(fm);
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

    bool is_mref_flipped(i64 i) const {
        return i >= static_cast<i32>(size() / 2);
    }

    bool is_mref_fwd(i64 i, bool is_rna) const {
        return is_mref_flipped(i) == is_rna;
    }

    i64 mref_to_pac(i64 mref) {
        if (is_mref_flipped(mref)) {
            return size() - mref + 1;
        }
        return mref;
    }

    i32 get_ref_id(i64 ref) {
        auto pac = mref_to_pac(ref);
        return bns_pos2rid(bns_, pac);
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


    RefCoord mrefs_to_ref_coord(i64 mref_start, i64 mref_end, bool read_fwd) const {

        bool flip = mref_start >= static_cast<i32>(size() / 2);

        i64 pac_st; //TODO rename
        if (flip) pac_st = size() - mref_end + 1;
        else pac_st = mref_start;

        i32 rid = bns_pos2rid(bns_, pac_st);
        if (rid < 0) return RefCoord();
        //assert(rid >= 0);

        auto name = std::string(bns_->anns[rid].name);
        auto start = pac_st - bns_->anns[rid].offset;
        auto end = start + (mref_end-mref_start);
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

    i64 ref_to_pac(std::string name, i64 coord) {
        i32 i;
        for (i = 0; i < bns_->n_seqs; i++) {
            if (strcmp(name.c_str(), bns_->anns[i].name) == 0)
                return bns_->anns[i].offset + coord;
        }
        return INT_MAX;
    } 

    i64 get_pac_shift(const std::string &ref_name) {
        auto rid = get_ref_id(ref_name);
        return bns_->anns[rid].offset;
    }

    bool pacseq_loaded() const {
        return pacseq_ != NULL;
    }

    std::pair<i64, i64> ref_to_mref(const std::string &ref_name, i64 st, i64 en, bool is_fwd, bool is_rna) {
        auto shift = get_pac_shift(ref_name);

        i64 pac_st = shift+st, 
            pac_en = shift+en;

        auto flip = is_fwd == is_rna;

        if (!flip) return {pac_st, pac_en};
        return {size() - pac_en, size() - pac_st};
    }

    i64 ref_to_mref(i32 rid, i64 ref_coord, bool is_fwd, bool is_rna) {
        auto shift = bns_->anns[rid].offset;
        auto pac_coord = shift+ref_coord;
        auto flip = is_fwd == is_rna;

        if (!flip) return pac_coord;
        return size() - pac_coord - 1;
    }

    i64 mref_to_ref(i64 mref) {
        i64 pac;
        if (is_mref_flipped(mref)) {
            pac = size() - mref - 1;
        } else {
            pac = mref;
        }
        i32 rid = bns_pos2rid(bns_, pac);
        return pac - bns_->anns[rid].offset;
    }

    u8 get_base(i64 pac, bool comp=false) {
        if (pac < 0 || pac > static_cast<i64>(size() / 2)) { //TODO better size definition
            throw std::out_of_range("Base out of range");
        }
        size_t i = pac >> 2,
               shift = ((pac & 3) ^ 3) << 1;
        return ((pacseq_[i] >> shift) & 3) ^ (comp | (comp<<1));
    }
    
    private:
    KmerType next_kmer(KmerType kmer, i64 pac, bool comp) {
        return ModelType::kmer_neighbor(kmer, get_base(pac, comp));
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

    std::vector<KmerType> get_kmers(i64 pac_start, i64 pac_end, bool rev, bool comp) {
        std::vector<KmerType> ret;
        if (!rev) {
            ret.push_back(get_kmer(pac_start, comp));
            for (auto i = pac_start+K; i < pac_end; i++) {
                next_kmer(ret, i, comp);
            }
        } else {
            ret.push_back(ModelType::kmer_rev(get_kmer(pac_end-K, comp)));
            for (auto i = pac_end-K-1; i >= pac_start; i--) {
                next_kmer(ret, i, comp);
            }
        }
        return ret;
    }

    std::vector<KmerType> get_kmers(const std::string &name, i64 start, i64 end, bool rev=false, bool comp=false) {
        i64 pac_start = ref_to_pac(name, start),
            pac_end = ref_to_pac(name, end);
        return get_kmers(pac_start, pac_end, rev, comp);
    }
    
    std::vector<KmerType> get_kmers(i64 mref_start, i64 mref_end, bool is_rna) {
        bool rev = is_mref_flipped(mref_end-1);
        bool comp = rev != is_rna;
        if (rev) {
            return get_kmers(size()-mref_end, size()-mref_start, true, comp);
        }
        return get_kmers(mref_start, mref_end, false, comp);
    }



    using FwdRevCoords = std::pair< std::vector<i64>, std::vector<i64> >;

    //Returns all FM index coordinates which translate into reference 
    //coordinates that overlap the specified range
    FwdRevCoords range_to_fms(std::string ref_name, i64 start, i64 end) {

        std::vector<i64> fwd_fms, rev_fms;

        auto ref_len = static_cast<i64>(size() / 2);

        auto slop = static_cast<int>( ceil(log(ref_len) / log(4)) );

        auto pac_min = ref_to_pac(ref_name, start),
             pac_max = pac_min + (end - start) - 1;

        i64 fwd_st;
        if (ref_len - pac_max > slop) {
            fwd_st = pac_max + slop;
        } else {
            fwd_st = ref_len - 1;
        }

        Range r = get_base_range(get_base(fwd_st));
        for (auto i = fwd_st-1; i >= pac_max && i <= fwd_st; i--) {
            r = get_neighbor(r, get_base(i));
        }

        for (auto f = r.start_; f <= r.end_; f++) {
            auto loc = fm_to_pac(f);
            if (loc == pac_max) {
                r = Range(f,f);
                break;
            }
        }

        fwd_fms.push_back(r.start_);
        for (auto i = pac_max-1; i >= pac_min && i < pac_max; i--) {
            r = get_neighbor(r, get_base(i));
            fwd_fms.push_back(r.start_);
        }

        i64 rev_st;
        if (pac_min > slop) {
            rev_st = pac_min - slop;
        } else {
            rev_st = 0;
        }

        r = get_base_range(BASE_COMP_B[get_base(rev_st)]);
        for (i64 i = rev_st+1; i <= pac_min; i++) {
            r = get_neighbor(r, BASE_COMP_B[get_base(i)]);
        }

        for (auto f = r.start_; f <= r.end_; f++) {
            auto loc = fm_to_mref(f);
            if (loc == pac_min) {
                r = Range(f,f);
                break;
            }
        }

        rev_fms.push_back(r.start_);
        for (auto i = pac_min+1; i <= pac_max; i++) {
            r = get_neighbor(r, BASE_COMP_B[get_base(i)]);
            rev_fms.push_back(r.start_);
        }

        return FwdRevCoords(rev_fms, fwd_fms);
    }

    class KmerSlice {
        public:
        KmerSlice(const u8 *pacseq, i64 st, i64 en) :
            pacseq_(pacseq),
            st_(st),
            en_(en),
            size_(en_-st_-K) {}

        KmerType operator[](i64 i) {
            i += st_;
            auto pst = i >> 2;
            u32 comb = *((u32 *) &pacseq_[pst]);
            auto shift = (i & 3) >> 1;
            return (KmerType) ( (comb >> shift) & ModelType::KMER_MASK );
        }

        i64 size() {
            return size_;
        }

        std::string to_str() {
            std::string str(size_, 'N');
            auto pst = st_ >> 2,
                pen = ((en_) >> 2)+1;
            u8 bst = (st_&3), ben;
            auto i = 0;
            for (auto j = pst; j < pen; j++) {
                ben = j == pen-1 ? (en_&3) : 4;
                for (u8 k = bst; k < ben; k++) {
                    str[i++] = BASE_CHARS[(pacseq_[j] >> ((k^3) << 1) ) & 3];
                }
                bst = 0;
            }
            return str;
        }

        private:
        const u8 *pacseq_;
        i32 st_, en_, pst_, pen_, bst_, ben_,
            size_;
    };

    //KmerSlice get_kmers

    #ifdef PYBIND

    #define PY_BWA_INDEX_METH(P) c.def(#P, &RefIndex<ModelType>::P);
    #define PY_BWA_INDEX_VEC(P) c.def(#P, py::vectorize(&RefIndex<ModelType>::P));

    static void pybind_defs(pybind11::module_ m, const std::string &suffix) {
        py::class_<RefIndex<ModelType>> c(m, ("RefIndex" + suffix).c_str());

        c.def(py::init<>());
        c.def(py::init<const std::string &>());
        c.def(py::init<const std::string &, bool>());
        c.def(py::init<const std::string &, bool, bool>());
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
        PY_BWA_INDEX_METH(fm_to_mref);
        PY_BWA_INDEX_METH(size);
        PY_BWA_INDEX_METH(mrefs_to_ref_coord);
        PY_BWA_INDEX_METH(get_seqs);
        PY_BWA_INDEX_METH(pacseq_loaded);
        PY_BWA_INDEX_METH(get_sa_loc);
        PY_BWA_INDEX_METH(get_ref_name);
        PY_BWA_INDEX_METH(get_ref_len);
        PY_BWA_INDEX_METH(get_pac_shift);
        PY_BWA_INDEX_METH(range_to_fms);
        PY_BWA_INDEX_METH(is_mref_fwd);
        PY_BWA_INDEX_METH(is_mref_flipped);
        PY_BWA_INDEX_VEC(get_base);
        c.def("get_ref_id", static_cast<i32 (RefIndex::*)(i64)> (&RefIndex::get_ref_id) );
        c.def("get_ref_id", static_cast<i32 (RefIndex::*)(const std::string &)> (&RefIndex::get_ref_id));
        //c.def("get_kmers_new", &RefIndex::get_kmers_new);
        c.def("ref_to_mref", static_cast<std::pair<i64, i64> (RefIndex::*)(const std::string &, i64, i64, bool, bool)> (&RefIndex::ref_to_mref));
        c.def("ref_to_mref", pybind11::vectorize(static_cast<i64 (RefIndex::*)(i32, i64, bool, bool)> (&RefIndex::ref_to_mref)));
        c.def("mref_to_ref", pybind11::vectorize(&RefIndex::mref_to_ref));
        //c.def("get_kmers", static_cast< std::vector<KmerType> (RefIndex::*)(i64, i64)> (&RefIndex::get_kmers) );
        
        c.def("get_kmers", 
            static_cast< std::vector<KmerType> (RefIndex::*)(i64, i64, bool, bool)> (&RefIndex::get_kmers),
            py::arg("pac_start"), py::arg("pac_end"), py::arg("rev"), py::arg("comp"));
        c.def("get_kmers", 
            static_cast< std::vector<KmerType> (RefIndex::*)(const std::string &, i64, i64, bool, bool)> (&RefIndex::get_kmers),
            py::arg("name"), py::arg("start"), py::arg("end"), py::arg("rev")=false, py::arg("comp")=false);
        c.def("get_kmers", 
            static_cast< std::vector<KmerType> (RefIndex::*)(i64, i64, bool)> (&RefIndex::get_kmers),
            py::arg("miref_start"), py::arg("miref_end"), py::arg("is_rna"));

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
