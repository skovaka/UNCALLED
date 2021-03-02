#ifndef _INCL_EVENT_PROFILER
#define _INCL_EVENT_PROFILER

#include <deque>
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

    std::vector<bool> get_full_mask(const std::vector<Event> &events) {
        std::vector<bool> mask;
        mask.reserve(events.size());

        reset();
        for (auto &e : events) {
            add_event(e);
            if (is_full()) {
                mask.push_back(to_mask_ == 0);
            }
        }
        
        while (mask.size() < events.size()) {
            if (to_mask_ == 0) {
                mask.push_back(true);
            } else {
                mask.push_back(false);
                to_mask_--;
            }
        }

        return mask;
    }

    #ifdef PYBIND

    #define PY_EVENT_PROFILER_METH(P) evpr.def(#P, &EventProfiler::P);
    #define PY_EVENT_PROFILER_PROP(P) evpr.def_property(#P, &EventProfiler::get_##P, &EventProfiler::set_##P);
    #define PY_EVENT_PROFILER_RPROP(P) evpr.def_property_readonly(#P, &EventProfiler::get_##P);
    #define PY_EVENT_PROFILER_PRM(P) prm.def_readonly(#P, &EventProfiler::Params::P);
    #define PY_ANNO_EVENT_VAL(P) ann.def_readonly(#P, &AnnoEvent::P);

    static void pybind_defs(pybind11::class_<EventProfiler> &evpr) {
        evpr.def(pybind11::init());
        PY_EVENT_PROFILER_METH(reset)
        PY_EVENT_PROFILER_METH(set_norm)
        PY_EVENT_PROFILER_METH(add_event)
        PY_EVENT_PROFILER_METH(is_full)
        PY_EVENT_PROFILER_METH(event_ready)
        PY_EVENT_PROFILER_METH(anno_event)
        PY_EVENT_PROFILER_METH(next_mean)
        PY_EVENT_PROFILER_METH(get_full_mask)

        pybind11::class_<AnnoEvent> ann(evpr, "AnnoEvent");
        PY_ANNO_EVENT_VAL(evt)
        PY_ANNO_EVENT_VAL(win_mean)
        PY_ANNO_EVENT_VAL(win_stdv)
        PY_ANNO_EVENT_VAL(mask)

        pybind11::class_<Params> prm(evpr, "Params");
        PY_EVENT_PROFILER_PRM(win_len)
        PY_EVENT_PROFILER_PRM(slop)
        PY_EVENT_PROFILER_PRM(win_stdv_min)
        PY_EVENT_PROFILER_PRM(win_stdv_range)
        PY_EVENT_PROFILER_PRM(win_mean_range)
    }
    #endif

};
#endif
