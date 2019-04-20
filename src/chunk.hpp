#ifndef CHUNK_TEST_HPP
#define CHUNK_TEST_HPP

#include <vector>
#include "util.hpp"

class Chunk {
    public:
    static void set_calibration(float digitisation, 
                                const std::vector<float> &offsets, 
                                const std::vector<float> &pa_ranges);

    Chunk();

    Chunk(const std::string &id, u16 channel, u32 number, u64 start_time, 
          const std::string &dtype, const std::string &raw_str);

    Chunk(const std::string &id, u16 channel, u32 number, u64 start_time, 
          const std::vector<float> &raw_data, u32 raw_st, u32 raw_len);
    Chunk(const Chunk &c);


    bool pop(std::vector<float> &raw_data);
    void swap(Chunk &c);
    void clear();

    float &operator[] (u32);

    bool empty() const;
    u64 get_start() const;
    u64 get_end() const;
    std::string get_id() const;
    u16 get_channel() const;
    u32 get_number() const;
    u32 size() const;
    bool is_calibrated() const;
    void print() const;


    private:
    float calibrate(i16 s);
    float calibrate(i32 s);

    //static float digitisation_;
    static std::vector<float> cal_offsets_, cal_coefs_;

    std::string id_;
    u16 channel_;
    u32 number_;
    u64 start_time_;
    std::vector<float> raw_data_;
    bool calibrated_;
    //std::vector<u32> chunk_classifications;
    //float median_before, median;

};

#endif
