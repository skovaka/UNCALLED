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
#include <unordered_map>
#include <iostream>
#include "event_detector.hpp"
#include "util.hpp"
#include "nt.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

typedef struct {
    KmerLen k;
    std::vector<float> means_stdvs;
} ModelPreset;

template<KmerLen KLEN>
class PoreModel {

    public:

    static const std::unordered_map<std::string, const std::vector<float> &> PRESETS;

    struct Params {
        std::string name;
        bool reverse, complement;
    };
    static const Params PRMS_DEF;

    Params PRMS;

    PoreModel(Params p) : PRMS(p) {
        loaded_ = false;

        kmer_count_ = kmer_count<KLEN>();

        kmer_means_.resize(kmer_count_);
        kmer_stdvs_.resize(kmer_count_);
        kmer_2vars_.resize(kmer_count_);
        lognorm_denoms_.resize(kmer_count_);

        if (p.name.empty()) return;

        auto vals = PRESETS.find(PRMS.name);
        if (vals != PRESETS.end()) {
            init_vals(vals->second);
        } else {
            init_tsv(PRMS.name);
        }
    }

    PoreModel() : PoreModel(PRMS_DEF) {
    }

    PoreModel(const std::vector<float> &means_stdvs, bool reverse, bool complement) 
        : PoreModel("", reverse, complement) {
        init_vals(means_stdvs);
    }

    PoreModel(const std::string &name, bool reverse=false, bool complement=false) : 
        PoreModel(Params({name, reverse, complement})) {}

    

    void init_vals(const std::vector<float> &vals) {
        model_mean_ = 0;

        u16 kmer = 0;
        for (u32 i = 0; i < vals.size(); i += 2) {
            float mean = vals[i],
                  stdv = vals[i+1];

            init_kmer(kmer, mean, stdv);
            
            kmer++;
            model_mean_ += mean;
        }

        model_mean_ /= kmer_count_;
        init_stdv();

        loaded_ = true;
    }
        
    void init_tsv(const std::string &model_fname) {
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

            init_kmer(kmer, lv_mean, lv_stdv);

            model_mean_ += lv_mean;
        }

        //Compute model level mean and stdv
        model_mean_ /= kmer_count_;
        init_stdv();

        loaded_ = true;
    }

    static bool is_preset(const std::string &name) {
        return PRESETS.find(name) != PRESETS.end();
    }

    static std::vector<std::string> get_preset_names() {
        std::vector<std::string> ret;
        for (auto p : PRESETS) {
            ret.push_back(p.first);
        }
        return ret;
    }

    float norm_pdf(float samp, u16 kmer) const {
        return (-pow(samp - kmer_means_[kmer], 2) / kmer_2vars_[kmer]) - lognorm_denoms_[kmer];
    }

    //TODO should be able to overload
    float match_prob_evt(const Event &evt, u16 kmer) const {
        return norm_pdf(evt.mean, kmer);
    }

    float abs_diff(float samp, u16 kmer) const {
        return std::abs(samp - kmer_means_[kmer]);
    }

    float model_mean() const {
        return model_mean_;
    }

    float model_stdv() const {
        return model_stdv_;
    }

    float kmer_current(u16 kmer) const {
        return kmer_means_[kmer];
    }

    float kmer_stdv(u16 kmer) const {
        return kmer_stdvs_[kmer];
    }

    bool is_loaded() const {
        return loaded_;
    }

    void calc_roc(std::vector<float> means, std::vector<u16> kmers, u32 n_threshs, bool prob_score) {

        float min_score = FLT_MAX,
              max_score = 0;
        std::vector<float> tp_scores(kmers.size());
        for (u64 i = 0; i < means.size(); i++) {
            if (prob_score) tp_scores[i] = -norm_pdf(means[i], kmers[i]);
            else tp_scores[i] = abs_diff(means[i], kmers[i]);
            min_score = std::min(min_score, tp_scores[i]);
            max_score = std::max(max_score, tp_scores[i]);
        }

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
                if (prob_score) score = -norm_pdf(means[i], k);
                else score = abs_diff(means[i], k);

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

    u32 get_kmer_count() const {
        return kmer_count_;
    }

    private:
    std::vector<float> kmer_means_, kmer_stdvs_, kmer_2vars_, lognorm_denoms_;
    float model_mean_, model_stdv_;
    u16 kmer_count_;
    bool loaded_, compl_;

    void init_stdv() {
        model_stdv_ = 0;

        for (u16 kmer = 0; kmer < kmer_count_; kmer++) {
            model_stdv_ += pow(kmer_means_[kmer] - model_mean_, 2);
        }

        model_stdv_ = sqrt(model_stdv_ / kmer_count_);
    }

    void init_kmer(u16 k, float mean, float stdv) {
        if (PRMS.reverse)  k = kmer_rev<KLEN>(k);
        if (PRMS.complement) k = kmer_comp<KLEN>(k);

        kmer_means_[k] = mean;
        kmer_stdvs_[k] = stdv;
        kmer_2vars_[k] = 2 * stdv * stdv;
        lognorm_denoms_[k] = log(sqrt(M_PI * kmer_2vars_[k]));
    }

    //void calc_roc_diff(std::vector<u16> kmers, std::vector<float> means, std::vector<float> threshs) {
    //    calc_roc(kmers, means, threshs, &PoreModel<KLEN>::abs_diff);
    //}

    #ifdef PYBIND

    #define PY_MODEL_DEF(P, D) c.def(#P, &PoreModel<KLEN>::P, D);
    #define PY_MODEL_DEFVEC(P, D) c.def(#P, py::vectorize(&PoreModel<KLEN>::P), D);
    #define PY_MODEL_PROP(P, D) c.def_property_readonly(#P, &PoreModel<KLEN>::P, D);
    #define PY_MODEL_PARAM(P, D) p.def_readwrite(#P, &PoreModel<KLEN>::Params::P, D);

    public:


    static void pybind_defs(py::module_ &m) {
        py::class_< PoreModel<KLEN> > c(m, "_PoreModel");

        py::class_<PoreModel<KLEN>::Params> p(c, "Params");
        p.def(py::init<>());
        p.def(py::init<PoreModel<KLEN>::Params>());
        PY_MODEL_PARAM(name, "Model preset name or TSV filename");
        PY_MODEL_PARAM(reverse, "Will reverse (flip) k-mer sequences if True");
        PY_MODEL_PARAM(complement, "Will complement k-mer sequences if True");

        c.def_readonly_static("PRMS_DEF", &PoreModel::PRMS_DEF);

        c.def(pybind11::init<const PoreModel<KLEN> &>());
        c.def(pybind11::init<PoreModel<KLEN>::Params>());
        c.def(pybind11::init<const std::string &, bool, bool>(), 
              py::arg("name"), py::arg("reverse")=PRMS_DEF.reverse, py::arg("complement")=PRMS_DEF.complement);
        c.def(pybind11::init<const std::vector<float> &, bool, bool>(), 
              py::arg("vals"), py::arg("reverse")=PRMS_DEF.reverse, py::arg("complement")=PRMS_DEF.complement);

        c.def_property_readonly("kmer_count", &PoreModel::get_kmer_count, "The number of k-mers in the model");
        PY_MODEL_PROP(model_mean, "The mean of all model k-mer currents");
        PY_MODEL_PROP(model_stdv, "The standard deviation of all model k-mer currents");

        c.def_static("is_preset", &PoreModel::is_preset, "List of model preset names");
        c.def_static("get_preset_names", &PoreModel::get_preset_names, "List of model preset names");

        c.def_property_readonly("means", 
            [](PoreModel<KLEN> &r) -> pybind11::array_t<float> {
                return pybind11::array_t<float>(r.kmer_means_.size(), r.kmer_means_.data());
        }, "The expected mean current of each k-mer");

        c.def_property_readonly("stdvs",
            [](PoreModel<KLEN> &r) -> pybind11::array_t<float> {
                return pybind11::array_t<float>(r.kmer_stdvs_.size(), r.kmer_stdvs_.data());
        }, "The expected standard devaition of each k-mer");

        PY_MODEL_DEFVEC(norm_pdf, "Returns the log probability that the current matches the k-mer based on the normal distibution probability density function");

        PY_MODEL_DEFVEC(abs_diff, "Returns the absolute difference between the observed and model current");

        c.def("__getitem__", py::vectorize(&PoreModel::kmer_current), "Alias for get_current()");

        c.def("__len__", &PoreModel::get_kmer_count, "Alias for kmer_count");

        PY_MODEL_DEF(calc_roc, "");

    }

    #endif
};

#endif
