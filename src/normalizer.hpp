#ifndef NORMALIZER_HPP
#define NORMALIZER_HPP

#include <vector>
#include "util.hpp"

class Normalizer {
    public:

    Normalizer();
    Normalizer(float tgt_mean, float tgt_stdv);

    void set_target(float tgt_mean, float tgt_stdv);

    void set_signal(const std::vector<float> &signal);

    float get_scale() const;
    float get_shift(float scale=0) const;

    float at(u32 i) const;

    float pop();
    bool push(float s);
    u32 skip_unread(u32 nkeep = 0);
    u32 unread_size() const;
    void reset(u32 buffer_size);

    bool empty() const;
    bool full() const;

    private:

    float tgt_mean_, tgt_stdv_;
    std::vector<float> signal_; //TODO: changed to float
    double mean_, varsum_;
    u32 n_, rd_, wr_;
    bool is_full_, is_empty_;
};

#endif
