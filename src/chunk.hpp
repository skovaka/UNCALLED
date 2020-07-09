#ifndef CHUNK_TEST_HPP
#define CHUNK_TEST_HPP

#include <vector>
#include "util.hpp"

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
    bool is_calibrated() const;
    void print() const;
    void set_start(u64 time);


    private:
    float calibrate(i16 s);
    float calibrate(i32 s);

    //static float digitisation_;
    static std::vector<float> cal_offsets_, cal_coefs_;

    std::string id_;
    u16 channel_idx_;
    u32 number_;
    u64 start_time_;
    std::vector<float> raw_data_;
    bool calibrated_;
    //std::vector<u32> chunk_classifications;
    //float median_before, median;

    friend bool operator< (const Chunk &r1, const Chunk &r2);
};

bool operator< (const Chunk &r1, const Chunk &r2);

#endif
