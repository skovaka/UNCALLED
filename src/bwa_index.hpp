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

#ifndef INCL_BWAFMI
#define INCL_BWAFMI

#include <string>
#include <climits>
#include <utility>
#include "util.hpp"
#include "bp.hpp"
#include "range.hpp"

#include "bwa/bwt.h"
#include "bwa/bntseq.h"
#include "bwa/utils.h"


template <KmerLen KLEN>
class BwaIndex {
    public:

    BwaIndex() :
        index_(NULL),
        bns_(NULL),
        pacseq_(NULL),
        klen_(KLEN),
        kmer_ranges_(kmer_count<KLEN>()),
        loaded_(false) {}

    BwaIndex(const std::string &prefix) {
        if (!prefix.empty()) load_index(prefix);
    }

    void load_index(const std::string &prefix) {
        std::string bwt_fname = prefix + ".bwt",
                    sa_fname = prefix + ".sa";

        index_ = bwt_restore_bwt(bwt_fname.c_str());
        bwt_restore_sa(sa_fname.c_str(), index_);
        bns_ = bns_restore(prefix.c_str());

        for (u16 k = 0; k < kmer_ranges_.size(); k++) {

            Range r = get_base_range(kmer_head<KLEN>(k));
            for (u8 i = 1; i < KLEN; i++) {
                r = get_neighbor(r, kmer_base<KLEN>(k, i));
            }

            kmer_ranges_[k] = r;
        }

        loaded_ = true;
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
    }

    Range get_neighbor(Range r1, u8 base) const {
        u64 os, oe;
        bwt_2occ(index_, r1.start_ - 1, r1.end_, base, &os, &oe);
        return Range(index_->L2[base] + os + 1, index_->L2[base] + oe);
    }

    Range get_kmer_range(u16 kmer) const {
        return kmer_ranges_[kmer];
    }

    Range get_base_range(u8 base) const {
        return Range(index_->L2[base], index_->L2[base+1]);
    }

    u64 sa(u64 i) const {
        return bwt_sa(index_, i);
    }

    u64 size() const {
        return index_->seq_len;
    }

    u64 translate_loc(u64 sa_loc, std::string &ref_name, u64 &ref_loc) const {
        i32 rid = bns_pos2rid(bns_, sa_loc);
        if (rid < 0) return 0;

        ref_name = std::string(bns_->anns[rid].name);
        ref_loc = sa_loc - bns_->anns[rid].offset;
        return bns_->anns[rid].len;
    }

    bool pacseq_loaded() const {
        return pacseq_ != NULL;

    }

    std::vector< std::pair<std::string, u64> > get_seqs() const {
        std::vector< std::pair<std::string, u64> > seqs;

        for (i32 i = 0; i < bns_->n_seqs; i++) {
            bntann1_t ann = bns_->anns[i];
            std::string name = std::string(ann.name);
            seqs.push_back( std::pair<std::string, u64>(name, ann.len) );
        }

        return seqs;
    }

    u8 get_base(u64 i) {
        return (pacseq_[i>>2] >> ( ((3^i)&3) << 1 )) & 3;
    }

    //void test() {
    //    BwaIndex<KmerLen::k5> i;
    //}

    //Range get_neighbor(Range range, u8 base) const;
    //Range get_base_range(u8 base) const;
    //Range get_kmer_range(u16 kmer) const;
    //u64 sa(u64 i) const;
    //u64 size() const;
    //u64 translate_loc(u64 sa_loc, std::string &ref_name, u64 &ref_loc) const;
    //std::vector< std::pair<std::string, u64> > get_seqs() const;
    //bool pacseq_loaded() const;
    //void load_pacseq();
    //u8 get_base(u64 i);

    private:
    bwt_t *index_;
    bntseq_t *bns_;
    u8 *pacseq_;
    KmerLen klen_;
    std::vector<Range> kmer_ranges_;
    bool loaded_;

};


#endif
