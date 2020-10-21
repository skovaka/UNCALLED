#ifndef _INCL_EVENT_PROFILER
#define _INCL_EVENT_PROFILER

#include "event_detector.hpp"
#include "normalizer.hpp"

typedef struct {
    Event evt;
    float win_mean;
    float win_stdv;
    u32 mask;
} AnnoEvent;

class EventProfiler {

    private:
    float norm_scale_{1}, 
          norm_shift_{0};

    Event next_evt_{0};
    float win_mean_, win_stdv_;

    std::deque<Event> events_;
    Normalizer window_;

    u32 total_count_{0};

    bool next_mask_{false}, is_full_{false};
    u32 to_mask_;
    const u32 WIN_MID;

    public: 
    typedef struct Params {
        u32 win_len, slop;
        float win_stdv_min;
        float win_stdv_range;
        float win_mean_range;
    } Params;

    static Params const PRMS_DEF;

    Params PRMS;

    std::vector<u32> mask_idx_map_;

    EventProfiler() : EventProfiler(PRMS_DEF) {};

    EventProfiler(Params p) : 
        WIN_MID(p.win_len / 2),
        PRMS(p) {
        window_.set_length(PRMS.win_len);
    }

    void reset() {
        window_.reset();
        events_.clear();
        next_evt_ = {0};
        is_full_ = false;
        to_mask_ = 0;
        
        //TODO only in debug mode?
        mask_idx_map_.clear();
        total_count_ = 0;
    }

    void set_norm(float scale, float shift) {
        norm_scale_ = scale;
        norm_shift_ = shift;
    }

    bool add_event(Event e) {
        window_.push(e.mean);
        events_.push_back(e);
        

        if (window_.unread_size() <= WIN_MID) return false;

        //float win_stdv = norm_scale_ * window_.get_stdv() + norm_shift_;
        win_mean_ = window_.get_mean();
        win_stdv_ = window_.get_stdv();
        if (win_stdv_ < PRMS.win_stdv_min) {
            to_mask_ = PRMS.win_len-1;
        } else if (to_mask_ > 0) {
            to_mask_--;
        }

        //TODO dynamic range bounds?

        if (window_.full()) {
            next_evt_ = events_.front();
            events_.pop_front();
            window_.pop();
            is_full_ = true;

            //TODO only in debug mode?
            if (to_mask_ == 0) {
                mask_idx_map_.push_back(total_count_);
            }
            total_count_ += 1;//TODO only in debug mode?
        }
        //window_.pop();

        return event_ready();
    }

    //TODO store midpoint
    //need to decide between "radius" or enforce odd
    bool is_full() {
        return is_full_;
    }

    bool event_ready() {
        return is_full_ && to_mask_ == 0;
    }

    AnnoEvent anno_event() {
        return {
            next_evt_, 
            win_mean_, 
            win_stdv_, 
            event_ready()
        };
    }

    float next_mean() {
        return next_evt_.mean;
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
