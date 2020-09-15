#ifndef _INCL_CHUNK
#define _INCL_CHUNK

#include <vector>
#include "util.hpp"

#ifdef PYBIND
#include "pybind11/pybind11.h"
#endif

class Chunk {
    public:
    Chunk();

    //Chunk(const Chunk &c);

    Chunk(const std::string &id, u16 channel, u32 number, u64 start_time, 
          const std::string &dtype, const std::string &raw_str);

    Chunk(const std::string &id, u16 channel, u32 number, u64 start_time, 
          const std::vector<float> &raw_data, u32 raw_st, u32 raw_len);

    bool pop(std::vector<float> &raw_data);
    void swap(Chunk &c);
    void clear();

    float &operator[] (u32);

    bool empty() const;
    u64 get_start() const;
    u64 get_end() const;
    std::string get_id() const;
    u16 get_channel() const;
    u16 get_channel_idx() const;
    u32 get_number() const;
    u32 get_raw_data() const;
    u32 size() const;
    void print() const;
    void set_start(u64 time);

    #ifdef PYBIND

    #define PY_CHUNK_METH(P) c.def(#P, &Chunk::P);
    #define PY_CHUNK_PROP(P) c.def_property(#P, &Chunk::get_##P, &Chunk::set_##P);
    #define PY_CHUNK_RPROP(P) c.def_property_readonly(#P, &Chunk::get_##P);

    static void pybind_defs(pybind11::class_<Chunk> &c) {
        c.def(pybind11::init<
            const std::string &, //id, 
            u16, u32, u64, //channel, number, start
            const std::string &, //dtype
            const std::string & //raw_str
        >());
        c.def(pybind11::init<
            const std::string &, //id, 
            u16, u32, u64, //channel, number, start
            const std::vector<float> &, //raw_data, 
            u32, u32 //raw_st, raw_len
        >());
        PY_CHUNK_METH(pop);
        PY_CHUNK_METH(swap);
        PY_CHUNK_METH(empty);
        PY_CHUNK_METH(print);
        PY_CHUNK_METH(size);
        PY_CHUNK_RPROP(channel);
        PY_CHUNK_RPROP(number);
    }

    #endif

    private:

    //static float digitisation_;
    static std::vector<float> cal_offsets_, cal_coefs_;

    std::string id_;
    u16 channel_idx_;
    u32 number_;
    u64 start_time_;
    std::vector<float> raw_data_;
    //std::vector<u32> chunk_classifications;
    //float median_before, median;

    friend bool operator< (const Chunk &r1, const Chunk &r2);
};

bool operator< (const Chunk &r1, const Chunk &r2);

#endif
