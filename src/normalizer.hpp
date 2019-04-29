#ifndef NORMALIZER_HPP
#define NORMALIZER_HPP

#include <vector>
#include "kmer_model.hpp"
#include "util.hpp"

class Normalizer {
    public:

    Normalizer();

    bool add_event(float newevt);
    float pop_event();
    NormParams get_params() const;
    u32 unread_size() const;
    u32 skip_unread(u32 nkeep = 0);
    void reset(u32 buffer_size);
    bool empty() const;
    bool full() const;

    //private:
    std::vector<double> events_;
    double mean_, varsum_;
    u32 n_, rd_, wr_;
    bool is_full_, is_empty_;
};

#endif
