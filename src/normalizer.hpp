#ifndef NORMALIZER_HPP
#define NORMALIZER_HPP

#include <vector>
#include "kmer_model.hpp"
#include "util.hpp"

class Normalizer {
    public:

    Normalizer(const KmerModel &model, 
               const EventParams &event_params, 
               u32 buffer_size);

    bool add_sample(float s);
    float pop_event();

    private:
    const KmerModel &model_;
    EventDetector detector_;
    std::vector<float> sum_, sumsq_;
    u32 n_, rd_, wr_;
};

#endif
