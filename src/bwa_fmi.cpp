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
#include "bwa_fmi.hpp"

BwaFMI::BwaFMI() {
    loaded_ = false;
}

BwaFMI::BwaFMI(const std::string &prefix) {
    std::string bwt_fname = prefix + ".bwt",
                sa_fname = prefix + ".sa";

	index_ = bwt_restore_bwt(bwt_fname.c_str());
	bwt_restore_sa(sa_fname.c_str(), index_);
    bns_ = bns_restore(prefix.c_str());

    loaded_ = true;
}

Range BwaFMI::get_neighbor(Range r1, u8 base) const {
    u64 os, oe;
    bwt_2occ(index_, r1.start_ - 1, r1.end_, base, &os, &oe);
    return Range(index_->L2[base] + os + 1, index_->L2[base] + oe);
}

Range BwaFMI::get_full_range(u8 base) const {
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
