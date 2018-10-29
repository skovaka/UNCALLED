#ifndef SELF_ALN_REF_HPP
#include <vector>
#include <string>
#include "util.hpp"

std::vector< std::vector<u64> > self_align(const std::string &bwa_prefix,
                                         const std::string fasta_fname,
                                         u32 sample_dist);

#endif
