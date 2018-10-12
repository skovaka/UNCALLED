#ifndef INCL_BWAFMI
#define INCL_BWAFMI

#include <string>
#include <climits>
#include "util.hpp"
#include "range.hpp"
#include "bwa/bwt.h"
#include "bwa/bntseq.h"

class BwaFMI {
    public:

    BwaFMI();

    BwaFMI(const std::string &prefix);

    Range get_neighbor(Range range, u8 base) const;

    Range get_full_range(u8 base) const;

    u64 sa(u64 i) const;

    u64 size() const;
    bntseq_t *bns_;

    private:
    bwt_t *index_;
    bool loaded_;
};

#endif
