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
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
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
#include "params.hpp"

Params PARAMS;

Params::Params() : mode(Mode::UNINIT){}

//Map constructor
void Params::init_map(
        const std::string &_bwa_prefix,
        const std::string &_model_fname,
        const std::string &_param_preset,
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
        float _min_top_conf) {
    PARAMS = 
        Params(Mode::MAP,_bwa_prefix,_model_fname,_param_preset,_seed_len,_min_aln_len,_min_rep_len,
         _max_rep_copy,_max_consec_stay,_max_paths,_max_events_proc,0,0,_evt_winlen1,
         _evt_winlen2,_threads,_num_channels,0,0,0,_evt_thresh1,_evt_thresh2,
         _evt_peak_height,_evt_min_mean,_evt_max_mean,_max_stay_frac,_min_seed_prob,
         _min_mean_conf,_min_top_conf,0,0,0,0,true,true);
}

void Params::init_realtime (
        const std::string &_bwa_prefix,
        const std::string &_model_fname,
        const std::string &_param_preset,
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
        float _sample_rate) {
    PARAMS =
       Params(Mode::REALTIME,_bwa_prefix,_model_fname,_param_preset,_seed_len,_min_aln_len,
       _min_rep_len,_max_rep_copy,_max_consec_stay,_max_paths,_max_events_proc,
       _max_chunks_proc,_evt_buffer_len,_evt_winlen1,_evt_winlen2,_threads,
       _num_channels,_chunk_len,_evt_batch_size,_evt_timeout,_evt_thresh1,
       _evt_thresh2,_evt_peak_height,_evt_min_mean,_evt_max_mean,_max_stay_frac,
       _min_seed_prob,_min_mean_conf,_min_top_conf,_max_chunk_wait,0,0,0,true,true);
    PARAMS.set_calibration(_offsets, _ranges, _digitisation);
    PARAMS.set_sample_rate(_sample_rate);
}

    //Simulate constructor
void Params::init_sim(
        const std::string &_bwa_prefix,
        const std::string &_model_fname,
        const std::string &_param_preset,
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
        bool  _sim_even,
        bool  _sim_odd) {
    PARAMS =
       Params(Mode::SIMULATE,_bwa_prefix,_model_fname,_param_preset,_seed_len,_min_aln_len,
       _min_rep_len,_max_rep_copy,_max_consec_stay,_max_paths,_max_events_proc,
       _max_chunks_proc,_evt_buffer_len,_evt_winlen1,_evt_winlen2,_threads,
       _num_channels,_chunk_len,_evt_batch_size,_evt_timeout,_evt_thresh1,
       _evt_thresh2,_evt_peak_height,_evt_min_mean,_evt_max_mean,_max_stay_frac,
       _min_seed_prob,_min_mean_conf,_min_top_conf,_max_chunk_wait,
       _sim_speed,_sim_st,_sim_en,_sim_even,_sim_odd);
}

Params::Params(Mode _mode,
               const std::string &_bwa_prefix,
               const std::string &_model_fname,
               const std::string &_param_preset,
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
               bool  _sim_even,
               bool  _sim_odd) :
    mode               (_mode),
    fmi                (_bwa_prefix),
    model              (_model_fname, true),
    event_params       ({_evt_winlen1,_evt_winlen2,
                         _evt_thresh1,_evt_thresh2,
                         _evt_peak_height,
                         _evt_min_mean,_evt_max_mean}),
    seed_len           (_seed_len),
    min_aln_len        (_min_aln_len),
    min_rep_len        (_min_rep_len),
    max_rep_copy       (_max_rep_copy),
    max_paths          (_max_paths),
    max_consec_stay    (_max_consec_stay),
    max_events_proc    (_max_events_proc),
    max_chunks_proc    (_max_chunks_proc),
    evt_buffer_len     (_evt_buffer_len),
    threads            (_threads),
    num_channels       (_num_channels),
    chunk_len          (_chunk_len),
    evt_batch_size     (_evt_batch_size),
    evt_timeout        (_evt_timeout),
    max_stay_frac      (_max_stay_frac),
    min_seed_prob      (_min_seed_prob),
    min_mean_conf      (_min_mean_conf),
    min_top_conf       (_min_top_conf),
    max_chunk_wait     (_max_chunk_wait),
    sim_speed          (_sim_speed),
    sim_st             (_sim_st),
    sim_en             (_sim_en),
    sim_even           (_sim_even),
    sim_odd            (_sim_odd),
    sample_rate        (4000),
    bp_per_sec         (450),
    calib_digitisation (0),
    calib_offsets      (_num_channels, 0), 
    calib_coefs        (_num_channels, 0){

    //TODO: exception handling
    std::ifstream param_file(_bwa_prefix + INDEX_SUFF);
    std::string param_line;

    char *param_preset = (char *) _param_preset.c_str();
    prob_threshes.resize(64);

    while (getline(param_file, param_line)) {
        char *param_name = strtok((char *) param_line.c_str(), "\t");
        char *fn_str = strtok(NULL, "\t");
        //char *path_str = strtok(NULL, "\t");
        if ( !_param_preset.empty() && strcmp(param_name, param_preset) ) {
                continue;
        }

        u8 fmbin = prob_threshes.size() - 1;
        char *prob_str;
        while ( (prob_str = strtok(fn_str, ",")) != NULL ) {
            fn_str = NULL;
            prob_threshes[fmbin] = atof(prob_str);
            fmbin--;
        }

        for (;fmbin < prob_threshes.size(); fmbin--) {
            prob_threshes[fmbin] = prob_threshes[fmbin+1];
        }

        //while ( (prob_str = strtok(path_str, ",")) != NULL ) {
        //    path_str = NULL;
        //    path_threshes.push_back(atof(prob_str));
        //}
    }

    u64 max_len = 0;
    kmer_fmranges = std::vector<Range>(model.kmer_count());
    for (u16 k = 0; k < model.kmer_count(); k++) {

        Range r = fmi.get_full_range(model.get_first_base(k));
        for (u8 i = 1; i < model.kmer_len(); i++) {
            r = fmi.get_neighbor(r, model.get_base(k, i));
        }

        kmer_fmranges[k] = r;

        if (r.length() > max_len) {
            max_len = r.length();
        }
    }


    bp_per_samp = bp_per_sec / sample_rate;
}


float Params::get_prob_thresh(u64 fmlen) const {
    return prob_threshes[__builtin_clzll(fmlen)];
}

float Params::get_path_thresh(u32 pathlen) const {
    return path_threshes[min(pathlen, path_threshes.size())-1];
}

float Params::get_source_prob() const {
    return prob_threshes.front();
}

bool Params::check_map_conf(u32 seed_len, float mean_len, float second_len) {
    return (min_mean_conf > 0 && seed_len / mean_len >= min_mean_conf) ||
           (min_top_conf > 0  && seed_len / second_len >= min_top_conf);
}

u16 Params::get_max_events(u16 event_i) const {
    if (event_i + evt_batch_size > max_events_proc) 
        return max_events_proc - event_i;
    return evt_batch_size;
}

void Params::set_sample_rate(float rate) {
    sample_rate = rate;
    bp_per_samp = bp_per_sec / sample_rate;
}

void Params::calibrate(u16 ch, std::vector<float> samples) {
    for (float &s : samples) s = calibrate(ch, s);
}

std::vector<float> Params::calibrate(u16 ch, std::vector<i16> samples) {
    std::vector<float> cal(samples.size());
    for (u32 i = 0; i < samples.size(); i++) {
        cal[i] = calibrate(ch, samples[i]);
    }
    return cal;
}

float Params::calibrate(u16 ch, float sample) {
    return (sample + calib_offsets[ch-1]) * calib_coefs[ch-1];
}

void Params::set_calibration(const std::vector<float> &offsets, 
                             const std::vector<float> &pa_ranges,
                             float digitisation) {
    if (digitisation > 0) calib_digitisation = digitisation;
    calib_offsets = offsets;
    calib_coefs.reserve(pa_ranges.size());
    for (float p : pa_ranges) calib_coefs.push_back(p / digitisation);
}

void Params::set_calibration(u16 channel,
                             float offsets, 
                             float pa_ranges,
                             float digitisation) {
    if (digitisation > 0) calib_digitisation = digitisation;
    calib_offsets[channel-1] = offsets;
    calib_coefs[channel-1] = pa_ranges / digitisation;
}
