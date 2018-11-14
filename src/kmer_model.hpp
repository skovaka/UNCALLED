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

#ifndef _INCL_KMER_MODEL
#define _INCL_KMER_MODEL

#include <array>
#include <utility>
#include "util.hpp"
#include "event_detector.hpp"

//TODO: move to normalizer
typedef struct NormParams {
    double shift, scale;
} NormParams;

class KmerModel {
    public:
    typedef std::vector<u16>::const_iterator neighbor_itr;

    KmerModel(std::string model_fname, bool complement);
    ~KmerModel();

    NormParams get_norm_params(const std::vector<Event> &events) const;
    void normalize(std::vector<Event> &events, NormParams norm={0, 0}) const;

    NormParams get_norm_params(const std::vector<float> &events) const;
    void normalize(std::vector<float> &raw, NormParams norm={0, 0}) const;

    u16 get_neighbor(u16 k, u8 i) const;
    std::pair<neighbor_itr, neighbor_itr> get_neighbors(u16 k_id) const {
        return std::pair<neighbor_itr, neighbor_itr>
                    (neighbors_[k_id].begin(), neighbors_[k_id].end());
    }

    bool event_valid(const Event &e) const;

    float event_match_prob(const Event &evt, u16 k_id) const;
    float event_match_prob(float evt, u16 k_id) const;

    float get_stay_prob(Event e1, Event e2) const; 

    inline u8 kmer_len() const {return k_;}
    inline u16 kmer_count() const {return kmer_count_;}

    u16 kmer_to_id(std::string kmer, u64 offset = 0) const;
    u16 kmer_comp(u16 kmer);

    u8 get_first_base(u16 k) const;
    u8 get_last_base(u16 k) const;
    u8 get_base(u16 kmer, u8 i) const;

    void parse_fasta (std::ifstream &fasta_in, 
                 std::vector<u16> &fwd_ids, 
                 std::vector<u16> &rev_ids) const;

    double *lv_means_, *lv_vars_x2_, *lognorm_denoms_, *sd_means_, *sd_stdvs_;

    double lambda_, model_mean_, model_stdv_;
    private:
    u8 k_;
    u16 kmer_count_;
    bool complement_;

    std::vector< std::vector<u16> > neighbors_;
    u16 *rev_comp_ids_;
};

#endif

