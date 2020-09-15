#ifndef _INCL_NORMALIZER
#define _INCL_NORMALIZER

#include <vector>
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#endif

class Normalizer {
    public:

    Normalizer();
    Normalizer(float tgt_mean, float tgt_stdv);

    void set_target(float tgt_mean, float tgt_stdv);
    void set_signal(const std::vector<float> &signal);

    float get_mean() const;
    float get_stdv() const;

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

    #ifdef PYBIND

    #define PY_NORM_METH(P) c.def(#P, &Normalizer::P);

    static void pybind_defs(pybind11::class_<Normalizer> &c) {
        c.def(pybind11::init());
        c.def(pybind11::init<float, float>());

        PY_NORM_METH(set_target);
        PY_NORM_METH(set_signal);
        PY_NORM_METH(get_mean);
        PY_NORM_METH(get_stdv);
        PY_NORM_METH(get_scale);
        PY_NORM_METH(get_shift);
        PY_NORM_METH(pop);
        PY_NORM_METH(push);
        PY_NORM_METH(skip_unread);
        PY_NORM_METH(unread_size);
        PY_NORM_METH(reset);
        PY_NORM_METH(empty);
        PY_NORM_METH(full);
    }

    #endif

    private:

    float tgt_mean_, tgt_stdv_;
    std::vector<float> signal_; //TODO: changed to float
    double mean_, varsum_;
    u32 n_, rd_, wr_;
    bool is_full_, is_empty_;
};

#endif
