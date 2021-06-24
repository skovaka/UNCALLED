#ifndef _INCL_NORMALIZER
#define _INCL_NORMALIZER

#include <vector>
#include "util.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

class Normalizer {
    public:
    struct Params {
        u32 len;
        float tgt_mean;
        float tgt_stdv;
    };

    static Params const PRMS_DEF;

    Params PRMS;

    Normalizer() : Normalizer(PRMS_DEF) {};
    Normalizer(Params p);
    Normalizer(float tgt_mean, float tgt_stdv);

    void set_length(u32 len);
    void set_target(float tgt_mean, float tgt_stdv);
    void set_signal(const std::vector<float> &signal);

    float kmer_current() const;
    float kmer_stdv() const;

    float get_scale() const;
    float get_shift() const;
    float get_shift(float scale) const;

    float at(u32 i) const;

    float pop();
    bool push(float s);
    u32 skip_unread(u32 nkeep = 0);
    u32 unread_size() const;
    void reset(u32 buffer_size = 0);

    bool empty() const;
    bool full() const;

    #ifdef PYBIND

    #define PY_NORM_METH(P, D) c.def(#P, &Normalizer::P, D);
    #define PY_NORM_PRM(P, D) p.def_readwrite(#P, &Normalizer::Params::P, D);

    public:

    static void pybind_defs(pybind11::class_<Normalizer> &c) {
        c.def(pybind11::init());
        c.def(pybind11::init<Params>());
        c.def(pybind11::init<float, float>());

        PY_NORM_METH(set_target, "Sets target mean and standard deviation");
        PY_NORM_METH(set_signal, "Sets full signal to normalize");
        PY_NORM_METH(set_length, "Sets the length of the rolling normalization window");
        PY_NORM_METH(kmer_current, "Returns the mean of the signal in the buffer");
        PY_NORM_METH(kmer_stdv, "Return the standard deviation of the signal in the buffer");
        PY_NORM_METH(get_scale, "Returns the scaling parameter required to normalize the signal in the buffer to the target mean/stdv");
        //PY_NORM_METH(get_shift, "Returns the shift parameter required to normalize the signal in the buffer to the target mean/stdv");

        c.def("get_shift", static_cast< float (Normalizer::*)() const> (&Normalizer::get_shift) );

        PY_NORM_METH(pop, "Returns the next unread normalized signal from the buffer");
        c.def("push", pybind11::vectorize(&Normalizer::push), "Adds new signal to the buffer");
        PY_NORM_METH(skip_unread, "Sets all signal in the buffer as read");
        PY_NORM_METH(unread_size, "Returns number of signals in the buffer which have not been read");
        PY_NORM_METH(reset, "Empties the buffer (buffer size and target mean/stdv are kept the same");
        PY_NORM_METH(empty, "Returns true if all signal has been read from the buffer");
        PY_NORM_METH(full, "Returns true if no more signal can be added without erasing the oldest signal");

        pybind11::class_<Params> p(c, "Params");
        PY_NORM_PRM(len, "The length of the normalization buffer")
        PY_NORM_PRM(tgt_mean, "Normalization target mean")
        PY_NORM_PRM(tgt_stdv, "Normalization target standard deviation")
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
