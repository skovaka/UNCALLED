#ifndef NORMALIZER_HPP
#define NORMALIZER_HPP

#include <vector>
#include "kmer_model.hpp"
#include "util.hpp"

class Normalizer {
    public:

    Normalizer(const KmerModel &model, 
               u32 buffer_size);

    bool add_event(float newevt);
    float pop_event();
    NormParams get_params() const;
    void reset(u32 buffer_size);
    bool empty() const;
    bool full() const;

    //private:
    const KmerModel &model_;
    std::vector<double> events_;
    double mean_, varsum_;
    u32 n_, rd_, wr_;
    bool is_full_, is_empty_;
};

#endif
