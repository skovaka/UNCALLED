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

#ifndef UNCALLED_OPTS_HPP
#define UNCALLED_OPTS_HPP

#include <iostream>
#include <vector>
#include "bwa_fmi.hpp"
#include "kmer_model.hpp"
#include "timer.hpp"

#define INDEX_SUFF ".uncl"

class UncalledOpts {
    public:

    enum Mode {MAP, SIMULATE};

    //Map constructor
    UncalledOpts(const std::string &bwa_prefix,
                 const std::string &model_fname,
                 u32 seed_len, 
                 u32 min_aln_len,
                 u32 min_rep_len, 
                 u32 max_rep_copy, 
                 u32 max_consec_stay,
                 u32 max_paths, 
                 u32 max_events_proc,
                 u32 evt_winlen1,
                 u32 evt_winlen2,
                 u16 threads,
                 float evt_thresh1,
                 float evt_thresh2,
                 float evt_peak_height,
                 float evt_min_mean,
                 float evt_max_mean,
                 float max_stay_frac,
                 float min_seed_prob, 
                 float min_mean_conf,
                 float min_top_conf);

    //Simulate constructor
    UncalledOpts(const std::string &bwa_prefix,
                 const std::string &model_fname,
                 u32 seed_len, 
                 u32 min_aln_len,
                 u32 min_rep_len, 
                 u32 max_rep_copy, 
                 u32 max_consec_stay,
                 u32 max_paths, 
                 u32 max_events_proc,
                 u32 max_chunks_proc,
                 u32 evt_buffer_len,
                 u32 evt_winlen1,
                 u32 evt_winlen2,
                 u16 threads,
                 u16 num_channels,
                 u16 chunk_len,
                 u16 evt_batch_size,
                 float evt_timeout,
                 float evt_thresh1,
                 float evt_thresh2,
                 float evt_peak_height,
                 float evt_min_mean,
                 float evt_max_mean,
                 float max_stay_frac,
                 float min_seed_prob, 
                 float min_mean_conf,
                 float min_top_conf,
                 float sim_speed);

    void init(const std::string &bwa_prefix,
                    const std::string &model_fname,
                    u32 seed_len, 
                    u32 min_aln_len,
                    u32 min_rep_len, 
                    u32 max_rep_copy, 
                    u32 max_consec_stay,
                    u32 max_paths, 
                    u32 max_events_proc,
                    u32 max_chunks_proc,
                    u32 evt_buffer_len,
                    u32 evt_winlen1,
                    u32 evt_winlen2,
                    u16 threads,
                    u16 num_channels,
                    u16 chunk_len,
                    u16 evt_batch_size,
                    float evt_timeout,
                    float evt_thresh1,
                    float evt_thresh2,
                    float evt_peak_height,
                    float evt_min_mean,
                    float evt_max_mean,
                    float max_stay_frac,
                    float min_seed_prob, 
                    float min_mean_conf,
                    float min_top_conf,
                    float sim_speed);

    
    u16 get_max_events(u16 event_i) const;
    float get_prob_thresh(u64 fm_length) const;
    float get_source_prob() const;

    Mode mode_;

    BwaFMI fmi_;
    KmerModel model_;
    EventParams event_params_;

    u32 seed_len_,
        min_rep_len_,
        max_rep_copy_,
        max_paths_,
        max_consec_stay_,
        min_aln_len_,
        max_events_proc_,
        max_chunks_proc_,
        evt_buffer_len_;

    u16 threads_,
        num_channels_,
        chunk_len_,
        evt_batch_size_;

    float evt_timeout_,
          max_stay_frac_,
          min_seed_prob_,
          min_mean_conf_,
          min_top_conf_,
          sim_speed_;

    std::vector<u64> evpr_lengths_;
    std::vector<float> evpr_threshes_;
    std::vector<Range> kmer_fmranges_;
};

#endif
