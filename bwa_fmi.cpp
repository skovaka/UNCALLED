#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <climits>
#include "bwa_fmi.hpp"

//void str_to_vec(const std::string &str, sdsl::int_vector<8> &vec) {
//    vec = sdsl::int_vector<8>(str.size());
//    for (size_t i = 0; i < str.size(); i++) {
//        vec[i] = char_to_base(str[i])+1;
//    }
//}

BwaFMI::BwaFMI() {
    loaded_ = false;
}

BwaFMI::BwaFMI(const std::string &prefix) {
    std::string bwt_fname = prefix + ".bwt",
                sa_fname = prefix + ".sa";

	index_ = bwt_restore_bwt(bwt_fname.c_str());
	bwt_restore_sa(sa_fname.c_str(), index_);

    loaded_ = true;
}

void BwaFMI::construct(const std::string &seq) {
    std::cerr << "Error: not implemented, use \"bwa index\"\n";
    loaded_ = false;
}

void BwaFMI::save(const std::string &filename) {
    std::cerr << "Error: not implemented, use \"bwa index\"\n";
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

