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

#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <climits>
#include <cmath>
#include "bwa_fmi.hpp"
#include "bwa/utils.h"


BwaFMI NULL_FMI;

BwaFMI::BwaFMI() {
    loaded_ = false;
    index_ = NULL;
    bns_ = NULL;
}

BwaFMI::BwaFMI(const std::string &prefix, const KmerModel &model) {
    pacseq_ = NULL;

    if (prefix.empty()) {
        loaded_ = false;
        index_ = NULL;
        bns_ = NULL;
    } else {
        std::string bwt_fname = prefix + ".bwt",
                    sa_fname = prefix + ".sa";

        index_ = bwt_restore_bwt(bwt_fname.c_str());
        bwt_restore_sa(sa_fname.c_str(), index_);
        bns_ = bns_restore(prefix.c_str());

        loaded_ = true;
    }


    //TESTING


	//u8 *pacseq;


    //for (u32 i = 0; i < 10; i++) {
    //    for (u8 j = 0; j < 4; j++) {
    //        u8 b = (pacseq[i] >> (j*2)) & 0x3;
    //        std::cout << (int) b << " ";
    //    }
    //    std::cout << "\n";
    //}

    //TESTING


    if (model.is_loaded()) {
        u64 max_len = 0;
        kmer_ranges_ = std::vector<Range>(model.kmer_count());
        for (u16 k = 0; k < model.kmer_count(); k++) {

            Range r = get_base_range(model.get_first_base(k));
            for (u8 i = 1; i < model.kmer_len(); i++) {
                r = get_neighbor(r, model.get_base(k, i));
            }

            kmer_ranges_[k] = r;

            if (r.length() > max_len) {
                max_len = r.length();
            }
        }
    }
}

void BwaFMI::destroy() {
    if (index_ != NULL) { 
        bwt_destroy(index_);
    }
    if (bns_ != NULL) { 
        bns_destroy(bns_);
    }
}

Range BwaFMI::get_neighbor(Range r1, u8 base) const {
    u64 os, oe;
    bwt_2occ(index_, r1.start_ - 1, r1.end_, base, &os, &oe);
    return Range(index_->L2[base] + os + 1, index_->L2[base] + oe);
}

Range BwaFMI::get_kmer_range(u16 kmer) const {
    return kmer_ranges_[kmer];
}

Range BwaFMI::get_base_range(u8 base) const {
    return Range(index_->L2[base], index_->L2[base+1]);
}

u64 BwaFMI::sa(u64 i) const {
    return bwt_sa(index_, i);
}

u64 BwaFMI::size() const {
    return index_->seq_len;
}

u64 BwaFMI::translate_loc(u64 sa_loc, std::string &ref_name, u64 &ref_loc) const {
    i32 rid = bns_pos2rid(bns_, sa_loc);
    if (rid < 0) return 0;

    ref_name = std::string(bns_->anns[rid].name);
    ref_loc = sa_loc - bns_->anns[rid].offset;
    return bns_->anns[rid].len;
}

bool BwaFMI::pacseq_loaded() const {
    return pacseq_ != NULL;

}

std::vector< std::pair<std::string, u64> > BwaFMI::get_seqs() const {
    std::vector< std::pair<std::string, u64> > seqs;

    for (i32 i = 0; i < bns_->n_seqs; i++) {
        bntann1_t ann = bns_->anns[i];
        std::string name = std::string(ann.name);
        seqs.push_back( std::pair<std::string, u64>(name, ann.len) );
    }

    return seqs;
}

void BwaFMI::load_pacseq() {
    if (!pacseq_loaded()) {
        //Copied from bwa/bwase.c
        pacseq_ = (u8*) calloc(bns_->l_pac/4+1, 1);
        err_fread_noeof(pacseq_, 1, bns_->l_pac/4+1, bns_->fp_pac);
    }   
}

u8 BwaFMI::get_base(u64 i) {
    return (pacseq_[i>>2] >> ( ((3^i)&3) << 1 )) & 3;
}


