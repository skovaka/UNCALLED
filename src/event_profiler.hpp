#ifndef _INCL_EVENT_PROFILER
#define _INCL_EVENT_PROFILER

class EventProfiler {

    public: 
    typedef struct Params {
        u32 n;
        float win_stdv_min;
        float win_stdv_range;
        float win_mean_range;
    } Params;

    static Params const PRMS_DEF;

    Params PRMS;

    EventProfiler() : Params(p) {};
    EventProfiler(Params p) :
        window_{

    void set_norm(float scale, float shift);

    void add_event(Event e) {
        if 
    }

    private:
    float norm_scale_{1}, 
          norm_shift_{0};
    Normalizer window_;
    std::deque<Event> win_evts_;

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
