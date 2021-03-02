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
#include <cmath>
#include <functional>
#include <cfloat>
#include <algorithm>
#include "event_detector.hpp"
#include "util.hpp"
#include "bp.hpp"

typedef struct {
    KmerLen k;
    std::vector<float> means_stdvs;
} ModelPreset;

template<KmerLen KLEN>
class PoreModel {

    private:
    std::vector<float> lv_means_, lv_stdvs_, lv_vars_x2_, lognorm_denoms_;
    float model_mean_, model_stdv_;
    u16 kmer_count_;
    bool loaded_, compl_;

    void init_stdv() {
        model_stdv_ = 0;

        for (u16 kmer = 0; kmer < kmer_count_; kmer++) {
            model_stdv_ += pow(lv_means_[kmer] - model_mean_, 2);
        }

        model_stdv_ = sqrt(model_stdv_ / kmer_count_);
    }

    void init_kmer(u16 k, float mean, float stdv) {
        lv_means_[k] = mean;
        lv_stdvs_[k] = stdv;
        lv_vars_x2_[k] = 2 * stdv * stdv;
        lognorm_denoms_[k] = log(sqrt(M_PI * lv_vars_x2_[k]));
    }

    public:

    PoreModel() 
        :  loaded_(false) {

        kmer_count_ = kmer_count<KLEN>();

        lv_means_.resize(kmer_count_);
        lv_stdvs_.resize(kmer_count_);
        lv_vars_x2_.resize(kmer_count_);
        lognorm_denoms_.resize(kmer_count_);
    }
    
    //PoreModel(const ModelPreset &p, bool cmpl) : PoreModel() {
    PoreModel(const std::vector<float> &means_stdvs, bool cmpl) 
        : PoreModel() {

        //if (p.k != KLEN) return;

        model_mean_ = 0;

        u16 kmer = 0;
        for (u32 i = 0; i < means_stdvs.size(); i += 2) {
            float mean = means_stdvs[i],
                  stdv = means_stdvs[i+1];

            if (cmpl) {
                init_kmer(kmer_comp<KLEN>(kmer), mean, stdv);
            } else { 
                init_kmer(kmer, mean, stdv);
            }
            
            kmer++;
            model_mean_ += mean;
        }

        model_mean_ /= kmer_count_;
        init_stdv();

        loaded_ = true;
    }

    //TODO: clean up IO
    //maybe load from toml and/or header file
    //make scripts for model to toml and header?
    PoreModel(const std::string &model_fname, bool cmpl) : PoreModel () {

        std::ifstream model_in(model_fname);

        if (!model_in.is_open()) {
            std::cerr << "Error: failed to open pore model file\n";
            throw std::runtime_error("Failed to open file");
        }

        std::string _;
        std::getline(model_in, _);

        //Variables for reading model
        std::string kmer_str, neighbor_kmer;
        u16 kmer;
        float lv_mean, lv_stdv;

        model_mean_ = 0;

        //Read and store rest of the model
        for (u32 i = 0; i < kmer_count_; i++) {
        //while (!model_in.eof()) {
        //
            if (model_in.eof()) {
                std::cerr << "Error: ran out of k-mers\n";
                return;
            }

            model_in >> kmer_str >> lv_mean >> lv_stdv; 

            //Get unique ID for the kmer
            kmer = str_to_kmer<KLEN>(kmer_str);

            if (kmer >= kmer_count_) {
                std::cerr << "Error: kmer '" << kmer << "' is invalid\n";
                return;
            }

            //Complement kmer if needed
            if (cmpl) {
                kmer = kmer_comp<KLEN>(kmer);
            }

            init_kmer(kmer, lv_mean, lv_stdv);

            model_mean_ += lv_mean;
        }

        //Compute model level mean and stdv
        model_mean_ /= kmer_count_;
        init_stdv();

        loaded_ = true;
    }

    float match_prob(float samp, u16 kmer) const {
        return (-pow(samp - lv_means_[kmer], 2) / lv_vars_x2_[kmer]) - lognorm_denoms_[kmer];
    }

    //TODO should be able to overload
    float match_prob_evt(const Event &evt, u16 kmer) const {
        return match_prob(evt.mean, kmer);
    }

    float match_diff(float samp, u16 kmer) const {
        return std::abs(samp - lv_means_[kmer]);
    }

    std::vector<float> match_probs(std::vector<float> means, std::vector<u16> kmers) const {
        std::vector<float> probs(means.size());
        for (u32 i = 0; i < means.size(); i++) {
            probs[i] = match_prob(means[i], kmers[i]);
        }
        return probs;
    }

    float get_means_mean() const {
        return model_mean_;
    }

    float get_means_stdv() const {
        return model_stdv_;
    }

    float get_mean(u16 kmer) const {
        return lv_means_[kmer];
    }

    float get_stdv(u16 kmer) const {
        return lv_stdvs_[kmer];
    }

    bool is_loaded() const {
        return loaded_;
    }

    void calc_roc(std::vector<float> means, std::vector<u16> kmers, u32 n_threshs, bool prob_score) {

        float min_score = FLT_MAX,
              max_score = 0;
        std::vector<float> tp_scores(kmers.size());
        for (u64 i = 0; i < means.size(); i++) {
            if (prob_score) tp_scores[i] = -match_prob(means[i], kmers[i]);
            else tp_scores[i] = match_diff(means[i], kmers[i]);
            min_score = std::min(min_score, tp_scores[i]);
            max_score = std::max(max_score, tp_scores[i]);
        }

        //
        //float min_score = tp_scores[0];
        //float max_score = tp_scores[static_cast<u32>(tp_scores.size()*0.99)];
        //
        std::sort(tp_scores.begin(), tp_scores.end());

        u64 dt = tp_scores.size() / (n_threshs-1);

        std::vector<float> threshs;
        for (u64 t = 0; t < tp_scores.size(); t += dt) {
            threshs.push_back(tp_scores[t]);
        }
        threshs.push_back(tp_scores.back());

        //float thresh = min_score, 
        //      dt = (max_score-min_score) / n_threshs;
        //for (u64 t = 0; t < n_threshs; t++) {
        //    threshs[t] = thresh;
        //    thresh += dt;
        //}

        std::vector<u64> tp_counts, fp_counts;
        tp_counts.resize(threshs.size());
        fp_counts.resize(threshs.size());

        u64 tpr_denom = 0,//means.size();
            fpr_denom = 0;

        for (u16 k = 0; k < kmer_count_; k++) {
            for (u64 i = 0; i < means.size(); i++) {
                bool is_true = kmers[i] == k;
                auto &tgt = is_true ? tp_counts : fp_counts;

                float score;
                if (prob_score) score = -match_prob(means[i], k);
                else score = match_diff(means[i], k);

                for (u64 t = 0; t < threshs.size(); t++) {
                    tgt[t] += score <= threshs[t];
                }

                if (is_true) tpr_denom += 1;
                else fpr_denom += 1;
            }
        }
        
        std::cout << "thresh\ttpr\tfpr\tmean_matches\n";
        for (u64 t = 0; t < threshs.size(); t++) {
            std::cout << threshs[t] << "\t"
                      << (static_cast<float>(tp_counts[t]) / tpr_denom) << "\t"
                      << (static_cast<float>(fp_counts[t]) / fpr_denom) << "\t"
                      << (static_cast<float>(fp_counts[t]+tp_counts[t]) / means.size()) << "\n";
        }
    }

    //void calc_roc_diff(std::vector<u16> kmers, std::vector<float> means, std::vector<float> threshs) {
    //    calc_roc(kmers, means, threshs, &PoreModel<KLEN>::match_diff);
    //}

    #ifdef PYBIND

    #define PY_PORE_MODEL_METH(P) c.def(#P, &PoreModel<KLEN>::P);

    static void pybind_defs(pybind11::class_<PoreModel<KLEN>> &c) {
        c.def(pybind11::init<const std::string &, bool>());
        PY_PORE_MODEL_METH(match_prob);
        PY_PORE_MODEL_METH(match_probs);
        PY_PORE_MODEL_METH(get_means_mean);
        PY_PORE_MODEL_METH(get_means_stdv);
        PY_PORE_MODEL_METH(get_mean);
        PY_PORE_MODEL_METH(get_stdv);
        PY_PORE_MODEL_METH(calc_roc);
    }

    #endif
};

#endif
