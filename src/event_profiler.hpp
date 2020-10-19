#ifndef _INCL_EVENT_PROFILER
#define _INCL_EVENT_PROFILER

#include "event_detector.hpp"
#include "normalizer.hpp"

class EventProfiler {

    private:
    float norm_scale_{1}, 
          norm_shift_{0};
    Normalizer window_;
    std::deque<Event> events_;

    public: 
    typedef struct Params {
        u32 win_len;
        float win_stdv_min;
        float win_stdv_range;
        float win_mean_range;
    } Params;

    static Params const PRMS_DEF;

    Params PRMS;

    EventProfiler() : PRMS(PRMS_DEF) {};

    EventProfiler(Params p) : PRMS(p) {
        window_.set_length(PRMS.win_len);
    }

    void reset() {
        window_.reset();
        events_.clear();
    }

    void set_norm(float scale, float shift) {
        norm_scale_ = scale;
        norm_shift_ = shift;
    }

    Event next_event(Event e) {
        window_.push(e.mean);
        events_.push_back(e);

        Event evt_out{0,0,0,0};

        //TODO store midpoint
        //need to decide between "radius" or enforce odd
        if (window_.unread_size() >= PRMS.win_len / 2) {

            //TODO reverse-normalize thresholds
            //float win_mean = norm_scale_ * window_.get_mean() + norm_shift_;
            //float win_stdv = norm_scale_ * window_.get_stdv() + norm_shift_;
            float win_stdv = window_.get_stdv();
            
            //TODO dynamic range bounds?
            if (win_stdv >= PRMS.win_stdv_min) {
                evt_out = events_.front();
            }
            events_.pop_front();
            window_.pop();
        }

        return evt_out;
    }

    Event pop_event() {
        Event e = events_.front();
        events_.pop_front();
        return e;
    }

    //#ifdef PYBIND

    //#define PY_CHUNK_METH(P) c.def(#P, &Chunk::P);
    //#define PY_CHUNK_PROP(P) c.def_property(#P, &Chunk::get_##P, &Chunk::set_##P);
    //#define PY_CHUNK_RPROP(P) c.def_property_readonly(#P, &Chunk::get_##P);

    //static void pybind_defs(pybind11::class_<Chunk> &c) {
    //    c.def(pybind11::init<
    //        const std::string &, //id, 
    //        u16, u32, u64, //channel, number, start
    //        const std::string &, //dtype
    //        const std::string & //raw_str
    //    >());
    //    c.def(pybind11::init<
    //        const std::string &, //id, 
    //        u16, u32, u64, //channel, number, start
    //        const std::vector<float> &, //raw_data, 
    //        u32, u32 //raw_st, raw_len
    //    >());
    //    PY_CHUNK_METH(pop);
    //    PY_CHUNK_METH(swap);
    //    PY_CHUNK_METH(empty);
    //    PY_CHUNK_METH(print);
    //    PY_CHUNK_METH(size);
    //    PY_CHUNK_RPROP(channel);
    //    PY_CHUNK_RPROP(number);
    //}

};

#endif
