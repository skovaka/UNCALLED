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
#include "util.hpp"
#include "nt.hpp"
#include "range.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

//From submods/bwa/bwtindex.c
#define BWA_BLOCK_SIZE 10000000

struct RefLoc {
    i32 ref_id;
    std::string ref_name;
    i64 ref_len, start, end;
    bool fwd;
};

//struct MirrorCoords {
//    u32
//};

template <KmerLen KLEN>
class BwaIndex {
    public:

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

    BwaIndex() :
        index_(NULL),
        bns_(NULL),
        pacseq_(NULL),
        klen_(KLEN),
        kmer_ranges_(kmer_count<KLEN>()),
        loaded_(false),
        size_(0) {}

    BwaIndex(const std::string &prefix, bool pacseq=false, bool bwt=true) : BwaIndex() {
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

        for (u16 k = 0; k < kmer_ranges_.size(); k++) {

            Range r = get_base_range(kmer_head<KLEN>(k));
            for (u8 i = 1; i < KLEN; i++) {
                r = get_neighbor(r, kmer_base<KLEN>(k, i));
            }

            kmer_ranges_[k] = r;
        }

        loaded_ = true;
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

    Range get_kmer_range(u16 kmer) const {
        return kmer_ranges_[kmer];
    }

    i64 get_kmer_count(u16 kmer) const {
        return kmer_ranges_[kmer].length();
    }

    Range get_base_range(u8 base) const {
        return Range(index_->L2[base], index_->L2[base+1]);
    }

    i64 fm_to_refmir(i64 fm) {
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

    bool is_pac_rev(i64 pac) const {
        return pac > static_cast<i32>(size() / 2);
    }

    bool pac_to_refmir(i64 pac) {
        if (is_pac_rev(pac)) {
            return size() - pac;
        }
        return pac;
    }

    i32 get_ref_id(i64 ref) {
        return bns_pos2rid(bns_, pac_to_refmir(ref));
    }

    std::string get_ref_name(u32 rid) {
        return bns_->anns[rid].name;
    }

    std::pair<u32, i32> get_ref_coord(i64 ref) {
        ref = pac_to_refmir(ref);

        auto rid = get_ref_id(ref);
        return {
            rid, 
            ref - bns_->anns[rid].offset
        };
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


    //auto ref_loc = fmi.translate_loc(seeds.ref_st_, seeds.ref_en_.end_ + KLEN, read_.PRMS.seq_fwd);
    RefLoc refmir_to_ref_bound(i64 sa_start, i64 sa_end, bool read_fwd) const {

        bool flip = sa_start >= static_cast<i32>(size() / 2);

        i64 pac_st; //TODO rename
        if (flip) pac_st = size() - sa_end + 1;
        else pac_st = sa_start;

        i32 rid = bns_pos2rid(bns_, pac_st);
        if (rid < 0) return {};
        //assert(rid >= 0);

        RefLoc ret {
            ref_id   : rid,
            ref_name : std::string(bns_->anns[rid].name),
            ref_len  : bns_->anns[rid].len,
            start    : pac_st - bns_->anns[rid].offset,
            end      : ret.start + (sa_end-sa_start),
            fwd      : (!flip && read_fwd) || (flip && !read_fwd)
        };

        return ret;
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

    std::pair<i64, i64> ref_to_refmir(const std::string &ref_name, i64 st, i64 en, bool is_fwd, bool is_rna) {
        auto shift = get_pac_shift(ref_name);

        i64 pac_st = shift+st, pac_en = shift+en;

        auto flip = is_fwd == is_rna;

        if (!flip) return {pac_st, pac_en};
        return {size() - pac_en, size() - pac_st};
    }

    i64 ref_to_refmir(i32 rid, i64 ref_coord, bool is_fwd, bool is_rna) {
        auto shift = bns_->anns[rid].offset;
        auto pac_coord = shift+ref_coord;
        auto flip = is_fwd == is_rna;

        if (!flip) return pac_coord;
        return size() - pac_coord;
    }

    i64 refmir_to_ref(i64 refmir) {
        i64 pac;
        if (is_pac_rev(refmir)) {
            pac = size() - refmir;
        } else {
            pac = refmir;
        }
        i32 rid = bns_pos2rid(bns_, pac);
        return pac - bns_->anns[rid].offset;
    }

    std::vector<u16> get_kmers(const std::string &nm, i64 st, i64 en) {
        i64 pac_st = ref_to_pac(nm, st),
            pac_en = ref_to_pac(nm, en);
        return seq_to_kmers<KLEN>(pacseq_, pac_st, pac_en);
    }

    std::vector<kmer_t> get_kmers(i64 mir_st, i64 mir_en, bool is_rna) {
        bool flip = is_pac_rev(mir_en);// >= size() / 2;

        bool fwd = flip == is_rna;

        auto st = mir_st, en = mir_en;
        if (flip) {
            st = size() - mir_en;
            en = size() - mir_st;
        }

        auto kmers = seq_to_kmers<KLEN>(pacseq_, st, en);

        std::vector<kmer_t> ret;
        ret.reserve(kmers.size());

        //TODO refactor
        if (flip) {
            for (auto itr = kmers.rbegin(); itr < kmers.rend(); itr++) {
                auto kmer = kmer_rev<KLEN>(*itr);
                if (!fwd) kmer = kmer_comp<KLEN>(kmer);
                ret.push_back(kmer);
            }
        } else {
            for (auto itr = kmers.begin(); itr < kmers.end(); itr++) {
                auto kmer = *itr;
                if (!fwd) kmer = kmer_comp<KLEN>(kmer);
                ret.push_back(kmer);
            }
        }

        return ret;
    }

    u8 get_base(u64 i) {
        return (pacseq_[i>>2] >> ( ((3^i)&3) << 1 )) & 3;
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
            auto loc = fm_to_refmir(f);
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
            size_(en_-st_-KLEN) {}

        u16 operator[](i64 i) {
            i += st_;
            auto pst = i >> 2;
            u32 comb = *((u32 *) &pacseq_[pst]);
            u8 shift = i & 3;
            return (u16) ( (comb >> ((16-KLEN)<<1)) & KMER_MASK );
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

    #define PY_BWA_INDEX_METH(P) c.def(#P, &BwaIndex<KLEN>::P);

    static void pybind_defs(pybind11::class_<BwaIndex<KLEN>> &c) {
        c.def(pybind11::init<>());
        c.def(pybind11::init<const std::string &>());
        c.def(pybind11::init<const std::string &, bool>());
        c.def(pybind11::init<const std::string &, bool, bool>());
        PY_BWA_INDEX_METH(create);
        PY_BWA_INDEX_METH(load_index);
        PY_BWA_INDEX_METH(bwt_loaded);
        PY_BWA_INDEX_METH(load_pacseq);
        PY_BWA_INDEX_METH(destroy);
        PY_BWA_INDEX_METH(get_neighbor);
        PY_BWA_INDEX_METH(get_kmer_range);
        PY_BWA_INDEX_METH(get_kmer_count);
        PY_BWA_INDEX_METH(get_base_range);
        PY_BWA_INDEX_METH(fm_to_pac);
        PY_BWA_INDEX_METH(fm_to_refmir);
        PY_BWA_INDEX_METH(size);
        PY_BWA_INDEX_METH(refmir_to_ref_bound);
        PY_BWA_INDEX_METH(get_seqs);
        PY_BWA_INDEX_METH(pacseq_loaded);
        PY_BWA_INDEX_METH(get_base);
        PY_BWA_INDEX_METH(get_sa_loc);
        PY_BWA_INDEX_METH(get_ref_coord);
        PY_BWA_INDEX_METH(get_ref_name);
        PY_BWA_INDEX_METH(get_ref_len);
        PY_BWA_INDEX_METH(get_pac_shift);
        PY_BWA_INDEX_METH(range_to_fms);
        c.def("get_ref_id", static_cast<i32 (BwaIndex::*)(i64)> (&BwaIndex::get_ref_id) );
        c.def("get_ref_id", static_cast<i32 (BwaIndex::*)(const std::string &)> (&BwaIndex::get_ref_id));
        //c.def("get_kmers_new", &BwaIndex::get_kmers_new);
        c.def("ref_to_refmir", static_cast<std::pair<i64,i64> (BwaIndex::*)(const std::string &, i64, i64, bool, bool)> (&BwaIndex::ref_to_refmir));
        c.def("ref_to_refmir", pybind11::vectorize(static_cast<i64 (BwaIndex::*)(i32, i64, bool, bool)> (&BwaIndex::ref_to_refmir)));
        c.def("refmir_to_ref", pybind11::vectorize(&BwaIndex::refmir_to_ref));
        //c.def("get_kmers", static_cast< std::vector<u16> (BwaIndex::*)(i64, i64)> (&BwaIndex::get_kmers) );
        c.def("get_kmers", static_cast< std::vector<u16> (BwaIndex::*)(const std::string &, i64, i64)> (&BwaIndex::get_kmers) );
        c.def("get_kmers", static_cast< std::vector<u16> (BwaIndex::*)(i64, i64, bool)> (&BwaIndex::get_kmers) );
    }

    #endif

    private:
    bwt_t *index_;
    bntseq_t *bns_;
    u8 *pacseq_;
    KmerLen klen_;
    std::vector<Range> kmer_ranges_;
    bool loaded_;
    size_t size_;
};



#endif
