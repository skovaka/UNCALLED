/* MIT License
 *
 * Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <iostream>
#include <vector>
#include "bwa_fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"

#define INDEX_SUFF ".uncl"

class Params {
    public:
    enum Mode {UNINIT, MAP, REALTIME, SIMULATE};

    Params();
    
    //Map constructor
    static void init_map (
        const std::string &_bwa_prefix,
        const std::string &_model_fname,
        u32 _seed_len, 
        u32 _min_aln_len,
        u32 _min_rep_len, 
        u32 _max_rep_copy, 
        u32 _max_consec_stay,
        u32 _max_paths, 
        u32 _max_events_proc,
        u32 _evt_winlen1,
        u32 _evt_winlen2,
        u16 _threads,
        u16 _num_channels,
        float _evt_thresh1,
        float _evt_thresh2,
        float _evt_peak_height,
        float _evt_min_mean,
        float _evt_max_mean,
        float _max_stay_frac,
        float _min_seed_prob, 
        float _min_mean_conf,
        float _min_top_conf);
    
    //Realtime constructor
    static void init_realtime (
        const std::string &_bwa_prefix,
        const std::string &_model_fname,
        u32 _seed_len, 
        u32 _min_aln_len,
        u32 _min_rep_len, 
        u32 _max_rep_copy, 
        u32 _max_consec_stay,
        u32 _max_paths, 
        u32 _max_events_proc,
        u32 _max_chunks_proc,
        u32 _evt_buffer_len,
        u32 _evt_winlen1,
        u32 _evt_winlen2,
        u16 _threads,
        u16 _num_channels,
        u16 _chunk_len,
        u16 _evt_batch_size,
        float _evt_timeout,
        float _evt_thresh1,
        float _evt_thresh2,
        float _evt_peak_height,
        float _evt_min_mean,
        float _evt_max_mean,
        float _max_stay_frac,
        float _min_seed_prob, 
        float _min_mean_conf,
        float _min_top_conf,
        float _max_chunk_wait,
        float _digitisation,
        std::vector<float> _offsets, 
        std::vector<float> _ranges,
        float sample_rate=4000);

    //Simulate constructor
    static void init_sim (
        const std::string &_bwa_prefix,
        const std::string &_model_fname,
        u32 _seed_len, 
        u32 _min_aln_len,
        u32 _min_rep_len, 
        u32 _max_rep_copy, 
        u32 _max_consec_stay,
        u32 _max_paths, 
        u32 _max_events_proc,
        u32 _max_chunks_proc,
        u32 _evt_buffer_len,
        u32 _evt_winlen1,
        u32 _evt_winlen2,
        u16 _threads,
        u16 _num_channels,
        u16 _chunk_len,
        u16 _evt_batch_size,
        float _evt_timeout,
        float _evt_thresh1,
        float _evt_thresh2,
        float _evt_peak_height,
        float _evt_min_mean,
        float _evt_max_mean,
        float _max_stay_frac,
        float _min_seed_prob, 
        float _min_mean_conf,
        float _min_top_conf,
        float _max_chunk_wait,
        float _sim_speed,
        float _sim_st,
        float _sim_en,
        bool  _sim_even);

    u16 get_max_events(u16 event_i) const;
    float get_prob_thresh(u64 fm_length) const;
    float get_source_prob() const;
    bool check_map_conf(u32 seed_len, float mean_len, float second_len);
    
    void set_sample_rate(float rate);
    void calibrate(u16 ch, std::vector<float> samples);
    std::vector<float> calibrate(u16 ch, std::vector<i16> samples);
    float calibrate(u16 ch, float sample);
    
    bool calibration_set(u16 channel = 0);

    void set_calibration(const std::vector<float> &offsets, 
                         const std::vector<float> &pa_ranges,
                         float digitisation = -1);

    void set_calibration(u16 channel,
                         float offsets, 
                         float pa_ranges,
                         float digitisation = -1);

    Mode mode;

    BwaFMI fmi;
    KmerModel model;
    EventParams event_params;

    u32 seed_len,
        min_aln_len,
        min_rep_len,
        max_rep_copy,
        max_paths,
        max_consec_stay,
        max_events_proc,
        max_chunks_proc,
        evt_buffer_len;

    u16 threads,
        num_channels,
        chunk_len,
        evt_batch_size;

    float evt_timeout,
          max_stay_frac,
          min_seed_prob,
          min_mean_conf,
          min_top_conf,
          max_chunk_wait,
          bp_per_samp,
          sim_speed,
          sim_st,
          sim_en;
    
    bool sim_even;


    std::vector<u64> evpr_lengths;
    std::vector<float> evpr_threshes;
    std::vector<Range> kmer_fmranges;

    float sample_rate, bp_per_sec;
    float calib_digitisation;
    std::vector<float> calib_offsets, calib_coefs;

    private: 

    Params(Mode _mode,
           const std::string &_bwa_prefix,
           const std::string &_model_fname,
           u32 _seed_len, 
           u32 _min_aln_len,
           u32 _min_rep_len, 
           u32 _max_rep_copy, 
           u32 _max_consec_stay,
           u32 _max_paths, 
           u32 _max_events_proc,
           u32 _max_chunks_proc,
           u32 _evt_buffer_len,
           u32 _evt_winlen1,
           u32 _evt_winlen2,
           u16 _threads,
           u16 _num_channels,
           u16 _chunk_len,
           u16 _evt_batch_size,
           float _evt_timeout,
           float _evt_thresh1,
           float _evt_thresh2,
           float _evt_peak_height,
           float _evt_min_mean,
           float _evt_max_mean,
           float _max_stay_frac,
           float _min_seed_prob, 
           float _min_mean_conf,
           float _min_top_conf,
           float _max_chunk_wait,
           float _sim_speed,
           float _sim_st,
           float _sim_en,
           bool  _sim_even);
};

extern Params PARAMS;

#endif
