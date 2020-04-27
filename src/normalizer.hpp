#ifndef NORMALIZER_HPP
#define NORMALIZER_HPP

#include <vector>
#include "pore_model.hpp"
#include "util.hpp"

class Normalizer {
    public:

    Normalizer();
    Normalizer(float tgt_mean, float tgt_stdv);

    void set_target(float tgt_mean, float tgt_stdv);

    void set_events(const std::vector<float> events);

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
