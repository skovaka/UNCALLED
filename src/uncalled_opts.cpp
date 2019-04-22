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

#include "pdqsort.h"
#include "uncalled_opts.hpp"

//Map constructor
UncalledOpts::UncalledOpts(
        const std::string &bwa_prefix,
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
        float min_top_conf) 
       : mode_(Mode::MAP),
         fmi_() {
    init(bwa_prefix,model_fname,seed_len,min_aln_len,min_rep_len,max_rep_copy,
         max_consec_stay,max_paths,max_events_proc,0,0,evt_winlen1,evt_winlen2,
         threads,0,0,0,0,evt_thresh1,evt_thresh2,evt_peak_height,evt_min_mean,
         evt_max_mean,max_stay_frac,min_seed_prob,min_mean_conf,min_top_conf,0);
}

    //Simulate constructor
UncalledOpts::UncalledOpts(
        const std::string &bwa_prefix,
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
        float sim_speed)
        : mode_(Mode::SIMULATE),
          fmi_() {
    init(bwa_prefix,model_fname,seed_len,min_aln_len,min_rep_len,max_rep_copy,
         max_consec_stay,max_paths,max_events_proc,max_chunks_proc,
         evt_buffer_len,evt_winlen1,evt_winlen2,threads,num_channels,chunk_len,
         evt_batch_size,evt_timeout,evt_thresh1,evt_thresh2,evt_peak_height,
         evt_min_mean,evt_max_mean,max_stay_frac,min_seed_prob,min_mean_conf,
         min_top_conf,sim_speed);
}

void UncalledOpts::init(const std::string &bwa_prefix,
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
                        float sim_speed) {

    fmi_             = BwaFMI(bwa_prefix);
    model_           = KmerModel(model_fname, true);
    seed_len_        = seed_len;
    min_aln_len_     = min_aln_len;
    min_rep_len_     = min_rep_len;
    max_rep_copy_    = max_rep_copy;
    max_consec_stay_ = max_consec_stay;
    max_paths_       = max_paths;
    max_events_proc_ = max_events_proc;
    max_chunks_proc_ = max_chunks_proc;
    evt_buffer_len_  = evt_buffer_len;
    threads_         = threads;
    num_channels_    = num_channels;
    chunk_len_       = chunk_len;
    evt_batch_size_  = evt_batch_size;
    evt_timeout_     = evt_timeout;
    max_stay_frac_   = max_stay_frac;
    min_seed_prob_   = min_seed_prob;
    min_mean_conf_   = min_mean_conf;
    min_top_conf_    = min_top_conf;
    sim_speed_       = sim_speed;
    event_params_    = {evt_winlen1,evt_winlen2,
                        evt_thresh1,evt_thresh2,
                        evt_peak_height,
                        evt_min_mean,evt_max_mean};
    
    //TODO: exception handling
    std::ifstream infile(bwa_prefix + INDEX_SUFF);
    float prob, frac;
    u64 fmlen = 0;
    infile >> prob >> frac;
    evpr_threshes_.push_back(prob);
    while (fmlen != 1) {
        infile >> fmlen >> prob >> frac;
        evpr_lengths_.push_back(fmlen);
        evpr_threshes_.push_back(prob);
    }

    kmer_fmranges_ = std::vector<Range>(model_.kmer_count());
    for (u16 k = 0; k < model_.kmer_count(); k++) {
        Range r = fmi_.get_full_range(model_.get_last_base(k));
        for (u8 i = model_.kmer_len()-2; i < model_.kmer_len(); i--) {
            r = fmi_.get_neighbor(r, model_.get_base(k, i));
        }
        kmer_fmranges_[k] = r;
    }
}

float UncalledOpts::get_prob_thresh(u64 fm_length) const {
    auto pr = evpr_threshes_.begin();
    for (auto len = evpr_lengths_.begin(); len != evpr_lengths_.end(); len++) {
        if (fm_length > *len) {
            break;
        }
        pr++;
    }
    return *pr;
}

float UncalledOpts::get_source_prob() const {
    return evpr_threshes_.front();
}

u16 UncalledOpts::get_max_events(u16 event_i) const {
    if (event_i + evt_batch_size_ > max_events_proc_) 
        return max_events_proc_ - event_i;
    return evt_batch_size_;
}

