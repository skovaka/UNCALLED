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

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<u16>);
PYBIND11_MAKE_OPAQUE(std::vector<u32>);
//PYBIND11_MAKE_OPAQUE(std::vector<u32>);

void pybind_pore_model_params(py::module_ &m);

#endif

using KmerLen = u8;

struct PoreModelParams {
    std::string name;
    bool reverse, complement;
    int k, shift;

};

struct PoreModelPreset {
    PoreModelParams prms;
    const std::vector<float> &vals;

    std::pair<std::string, PoreModelPreset> map() const {
        return {prms.name, *this};
    }
};

extern const std::unordered_map<std::string, PoreModelPreset> PORE_MODEL_PRESETS;

extern const PoreModelParams PORE_MODEL_PRMS_DEF;

using PresetMap = std::unordered_map<std::string, const std::vector<float> &>;

template<KmerLen K, typename KmerType=typename std::conditional<(K < 8), u16, u32>::type>
class PoreModel {

    public:
    static constexpr KmerType KMER_MASK = (1 << (2*K)) - 1,
                     KMER_COUNT = static_cast<KmerType>(pow(BASE_COUNT, K));
    static constexpr KmerLen KMER_LEN = K;
    //static const KmerLen SHIFT;

    using kmer_t = KmerType;

    static std::unordered_map<std::string, PoreModelPreset &> MODEL_PRESETS;

    //static bool PRESETS_LOADED;
    //static const PresetMap PRESETS;

    //inline static bool load_presets() {
    //    if (!PRESETS_LOADED) {
    //        for (auto &preset : PORE_MODEL_PRESETS) {
    //            if (preset.prms.k == K) {
    //                MODEL_PRESETS[preset.prms.name] = preset;
    //            }
    //        }
    //        return (PRESETS_LOADED = true);
    //    }
    //    return false;
    //}

    //static bool is_preset(const std::string &name) {
    //    return PRESETS.find(name) != PRESETS.end();
    //}

    //static std::vector<std::string> get_preset_names() {
    //    std::vector<std::string> ret;
    //    for (auto p : PRESETS) {
    //        ret.push_back(p.first);
    //    }
    //    return ret;
    //}

    PoreModelParams PRMS;

    std::vector<float> kmer_means_, kmer_stdvs_, kmer_2vars_, lognorm_denoms_;
    float model_mean_, model_stdv_;
    bool loaded_, compl_;

    PoreModel(PoreModelParams p) : PRMS(p) {
        loaded_ = false;

        kmer_means_.resize(KMER_COUNT);
        kmer_stdvs_.resize(KMER_COUNT);
        kmer_2vars_.resize(KMER_COUNT);
        lognorm_denoms_.resize(KMER_COUNT);


        if (p.name.empty()) return;

        //load_presets();
        auto preset = PORE_MODEL_PRESETS.find(PRMS.name);

        if (preset != PORE_MODEL_PRESETS.end()) {
            init_vals(preset->second.vals);
        } else {
            init_tsv(PRMS.name);
        }
    }

    PoreModel() : PoreModel(PORE_MODEL_PRMS_DEF) {
    }

    PoreModel(const std::vector<float> &means_stdvs, bool reverse, bool complement) 
        : PoreModel("", reverse, complement) {
        init_vals(means_stdvs);
    }

    PoreModel(const std::string &name, bool reverse=false, bool complement=false) : 
            PoreModel(PoreModelParams({name, reverse, complement})) {

    }

    void init_vals(const std::vector<float> &vals) {
        model_mean_ = 0;

        KmerType kmer = 0;
        for (u32 i = 0; i < vals.size(); i += 2) {
            float mean = vals[i],
                  stdv = vals[i+1];

            init_kmer(kmer, mean, stdv);
            
            kmer++;
            model_mean_ += mean;
        }

        model_mean_ /= KMER_COUNT;
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
        KmerType kmer;
        float lv_mean, lv_stdv;

        model_mean_ = 0;

        //Read and store rest of the model
        for (u32 i = 0; i < KMER_COUNT; i++) {
            if (model_in.eof()) {
                std::cerr << "Error: ran out of k-mers\n";
                return;
            }

            model_in >> kmer_str >> lv_mean >> lv_stdv; 

            //Get unique ID for the kmer
            kmer = str_to_kmer(kmer_str);

            if (kmer >= KMER_COUNT) {
                std::cerr << "Error: kmer '" << kmer << "' is invalid\n";
                return;
            }

            init_kmer(kmer, lv_mean, lv_stdv);

            model_mean_ += lv_mean;
        }

        //Compute model level mean and stdv
        model_mean_ /= KMER_COUNT;
        init_stdv();

        loaded_ = true;
    }

    float norm_pdf(float samp, KmerType kmer) const {
        return (-pow(samp - kmer_means_[kmer], 2) / kmer_2vars_[kmer]) - lognorm_denoms_[kmer];
    }

    //TODO should be able to overload
    float match_prob_evt(const Event &evt, KmerType kmer) const {
        return norm_pdf(evt.mean, kmer);
    }

    float abs_diff(float samp, KmerType kmer) const {
        return std::abs(samp - kmer_means_[kmer]);
    }

    float z_score(float samp, KmerType kmer) const {
        return std::abs(samp - kmer_means_[kmer]) / kmer_stdvs_[kmer];
    }

    float model_mean() const {
        return model_mean_;
    }

    float model_stdv() const {
        return model_stdv_;
    }

    float kmer_current(KmerType kmer) const {
        return kmer_means_[kmer];
    }

    float kmer_stdv(KmerType kmer) const {
        return kmer_stdvs_[kmer];
    }

    bool is_loaded() const {
        return loaded_;
    }

    void calc_roc(std::vector<float> means, std::vector<KmerType> kmers, u32 n_threshs, bool prob_score) const {

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

        for (KmerType k = 0; k < KMER_COUNT; k++) {
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
        return KMER_COUNT;
    }

    void init_stdv() {
        model_stdv_ = 0;

        for (KmerType kmer = 0; kmer < KMER_COUNT; kmer++) {
            model_stdv_ += pow(kmer_means_[kmer] - model_mean_, 2);
        }

        model_stdv_ = sqrt(model_stdv_ / KMER_COUNT);
    }

    void init_kmer(KmerType k, float mean, float stdv) {
        if (PRMS.reverse)  k = kmer_rev(k);
        if (PRMS.complement) k = kmer_comp(k);

        kmer_means_[k] = mean;
        kmer_stdvs_[k] = stdv;
        kmer_2vars_[k] = 2 * stdv * stdv;
        lognorm_denoms_[k] = log(sqrt(M_PI * kmer_2vars_[k]));
    }

    static char base_to_char(u8 base) {
        return BASE_CHARS[base];
    }

    static u8 base_comp(u8 base) {
        return BASE_COMP_B[base];
    }

    static KmerType str_to_kmer(const std::string &kmer, u32 offs=0) {
        KmerType index = BASE_BYTES[(u8) kmer[offs]];
        for (u8 i = 1; i < (u8) K; i++) {
            index = (index << 2) | BASE_BYTES[(u8) kmer[offs+i]];
        }
        return index;
    }

    static KmerType kmer_comp(KmerType kmer) {
        return kmer ^ KMER_MASK;
    }

    private:
    static u16 _kmer_rev(u16 kmer) {
        kmer = ( (kmer >> 2 & 0x3333) | (kmer & 0x3333) << 2 );
        kmer = ( (kmer >> 4 & 0x0F0F) | (kmer & 0x0F0F) << 4 );
        kmer = ( (kmer >> 8 & 0x00FF) | (kmer & 0x00FF) << 8 );
        return kmer >> (2 * (8 - K));
    }

    static u32 _kmer_rev(u32 kmer) {
        kmer = ( (kmer >> 2  & 0x33333333) | (kmer & 0x33333333) << 2 );
        kmer = ( (kmer >> 4  & 0x0F0F0F0F) | (kmer & 0x0F0F0F0F) << 4 );
        kmer = ( (kmer >> 8  & 0x00FF00FF) | (kmer & 0x00FF00FF) << 8 );
        kmer = ( (kmer >> 16 & 0x0000FFFF) | (kmer & 0x0000FFFF) << 16 );
        return kmer >> (2 * (16 - K));
    }

    public:
    static KmerType kmer_rev(KmerType kmer) {
        return _kmer_rev(kmer);
    }

    static KmerType kmer_revcomp(KmerType kmer) {
        return kmer_rev(~kmer);
    }

    static std::vector<KmerType> kmers_revcomp(const std::vector<KmerType> &kmers) {
        std::vector<KmerType> rev;
        rev.reserve(kmers.size());
        for (auto k = kmers.rbegin(); k != kmers.rend(); k++) {
            rev.push_back(kmer_revcomp(*k));
        }
        return rev;
    }

    static u8 kmer_head(KmerType kmer) {
        return (u8) ((kmer >> (2*( (u8) K ) - 2)) & 0x3);
    }

    static KmerType kmer_neighbor(KmerType kmer, u8 i) {
        return ((kmer << 2) & KMER_MASK) | i; 
    }

    static KmerType set_kmer_base(KmerType kmer, u8 i, u8 b) {
        auto shift = 2*i;
        auto mask = ~(0x3 << shift);
        return (kmer & mask) | (b << shift); 
    }

    static u8 kmer_base(KmerType kmer, u8 i) {
        return (u8) ((kmer >> (2 * ((KmerType)K-i-1))) & 0x3);
    }

    static u8 kmer_base_count(KmerType kmer, u8 base) {
        u8 count = 0;
        for (auto i = 0; i < K; i++) {
            count += kmer_base(kmer, i) == base;
        }
        return count;
    }

    static std::string kmer_to_str(KmerType kmer) {
        std::string s(K, 'N');
        for (u8 i = 0; i < K; i++) {
            s[i] = BASE_CHARS[kmer_base(kmer, i)];
        }
        return s;
    }

    static ValArray<KmerType> pacseq_to_kmers(u8 *seq, u64 st, u64 en) {
        ValArray<KmerType> ret(en - st - K + 1);

        u64 pst = st >> 2,
            pen = ((en) >> 2)+1;

        u64 i = 0;
        KmerType kmer = 0;
        u8 bst = (st&3), ben;

        //TODO could optimize by splitting into init and extend loops
        for (u64 j = pst; j < pen; j++) {
            ben = j == pen-1 ? (en&3) : 4;
            for (u8 k = bst; k < ben; k++) {
                kmer = kmer_neighbor(kmer, (seq[j] >> ((K^3) << 1) ) & 3);
                if (i >= K-1) {
                    ret[i] = kmer;
                }
                i++;
            }
            bst = 0;
        }

        return ret;
    }

    #ifdef PYBIND

    #define PY_MODEL_DEF(P, D) c.def(#P, &PoreModel<K,KmerType>::P, D);
    #define PY_MODEL_DEF_STATIC(P, D) c.def_readwrite_static(#P, &PoreModel<K,KmerType>::P, D);
    #define PY_MODEL_DEF_CONSTANT(P, D) c.def_readonly_static(#P, &PoreModel<K,KmerType>::P, D);
    #define PY_MODEL_DEFVEC(P, D) c.def(#P, py::vectorize(&PoreModel<K,KmerType>::P), D);
    #define PY_MODEL_PROP(P, D) c.def_property_readonly(#P, &PoreModel<K,KmerType>::P, D);
    #define PY_MODEL_PARAM(P, D) p.def_readwrite(#P, &PoreModelParams::P, D);

    public:

	using KmerTypePy = std::array<char, K+1>;

	static KmerTypePy kmer_to_arr(KmerType kmer) {
		KmerTypePy ret{0};
		for (size_t i = 0; i < K; i++) {
			ret[i] = BASE_CHARS[kmer_base(kmer, i)];
		}
		return ret;
    }

    static KmerType str_to_kmer(const KmerTypePy &kmer, u32 offs) {
		auto s = std::string(kmer.data());
		return str_to_kmer(s, offs);
	}



    static void pybind_defs(py::module_ &m, const std::string &suffix) {
        using Class = PoreModel<K,KmerType>;

        py::class_<Class> c(m, ("PoreModel" + suffix).c_str());

        c.def(pybind11::init<const Class &>());
        c.def(pybind11::init<PoreModelParams>());
        c.def(pybind11::init<const std::string &, bool, bool>(), 
              py::arg("name"), py::arg("reverse")=PORE_MODEL_PRMS_DEF.reverse, py::arg("complement")=PORE_MODEL_PRMS_DEF.complement);
        c.def(pybind11::init<const std::vector<float> &, bool, bool>(), 
              py::arg("vals"), py::arg("reverse")=PORE_MODEL_PRMS_DEF.reverse, py::arg("complement")=PORE_MODEL_PRMS_DEF.complement);
        
        //c.attr("KmerArray") = m.attr("PyArray;

        c.def_property_readonly("kmer_count", &Class::get_kmer_count, "The number of k-mers in the model");
        c.def_readwrite("PRMS", &Class::PRMS, "Class parameters");
        PY_MODEL_PROP(model_mean, "The mean of all model k-mer currents");
        PY_MODEL_PROP(model_stdv, "The standard deviation of all model k-mer currents");

        //c.def_static("is_preset", &Class::is_preset, "List of model preset names");
        //c.def_static("get_preset_names", &Class::get_preset_names, "List of model preset names");

        c.def_property_readonly("means", 
            [](Class &r) -> pybind11::array_t<float> {
                return pybind11::array_t<float>(r.kmer_means_.size(), r.kmer_means_.data());
        }, "The expected mean current of each k-mer");

        c.def_property_readonly("stdvs",
            [](Class &r) -> pybind11::array_t<float> {
                return pybind11::array_t<float>(r.kmer_stdvs_.size(), r.kmer_stdvs_.data());
        }, "The expected standard devaition of each k-mer");

        PY_MODEL_DEFVEC(norm_pdf, "Returns the log probability that the current matches the k-mer based on the normal distibution probability density function");

        PY_MODEL_DEFVEC(abs_diff, "Returns the absolute difference between the observed and model current");

        c.def("__getitem__", py::vectorize(&Class::kmer_current), "Alias for get_current()");

        c.def("__len__", &Class::get_kmer_count, "Alias for kmer_count");

        //c.def_readonly_static("K", &Class::KMER_LEN);
        //c.def_readonly_static("KMER_COUNT", &Class::KMER_COUNT);
        c.attr("K") = pybind11::cast(Class::KMER_LEN);
        c.attr("KMER_COUNT") = pybind11::cast(Class::KMER_COUNT);

        c.def_static("str_to_kmer",
                py::vectorize(static_cast< KmerType (*) (const KmerTypePy &, u32)>(&Class::str_to_kmer)), 
                py::arg("kmer"), py::arg("offs")=0);

        c.def_static("kmer_to_str",  &Class::kmer_to_str, "Convert binary k-mers to strings");
        c.def_static("kmer_to_arr",  py::vectorize(&Class::kmer_to_arr));

        c.def_static("base_to_char", py::vectorize(&Class::base_to_char), "Convert base index to base character"); 
        c.def_static("base_comp", py::vectorize(&Class::base_comp), "Returns the complement of a base index");
        
        c.def_static("kmer_rev",      py::vectorize(&Class::kmer_rev), "Reverses binary k-mers (not complement)");
        c.def_static("kmer_comp",     py::vectorize(&Class::kmer_comp), "Complements binary k-mers (not reverse)");
        c.def_static("kmer_revcomp",  py::vectorize(&Class::kmer_revcomp), "Reverse complements binary k-mers");
        c.def_static("kmer_head",     py::vectorize(&Class::kmer_head), "Returns the first base index in binary k-mers");
        c.def_static("kmer_base",     py::vectorize(&Class::kmer_base), "Returns the base at the specified index in a binary k-mer");
        c.def_static("set_kmer_base",     py::vectorize(&Class::set_kmer_base), "Sets a binary k-mer base");
        c.def_static("kmer_base_count",     py::vectorize(&Class::kmer_base_count), "Returns the number of occurances of the specified base index in the binary k-mer");
        c.def_static("kmer_neighbor", py::vectorize(&Class::kmer_neighbor), "Returns the binary k-mer shifted left with the specified base appended");


    }

    #endif
};

template <typename ModelType>
struct Sequence {//: public DataFrame<typename ModelType::kmer_t, float, u8> {
    using KmerType = typename ModelType::kmer_t;
    //using super = DataFrame<KmerType, u8, float>;

    static const KmerType K = ModelType::KMER_LEN;

    const ModelType &model;
    i32 ref_id;
    IntervalIndex<i64> coords;
    bool is_fwd; //TODO infer from mref coords

    //static constexpr typename super::NameArray names = {"ref", "start", "end"}; 
    //typename super::template ColType<0> &kmer = std::get<0>(super::data_);   
    //typename super::template ColType<1> &current = std::get<1>(super::data_);   
    //typename super::template ColType<2> &base = std::get<2>(super::data_);   
    //ValArray<u8> base;

    ValArray<KmerType> kmer; 
    ValArray<float> current;   

    Sequence(const ModelType &model_, size_t length) : 
        model(model_), 
        ref_id(-1),
        coords({{0,static_cast<i64>(length)}}), 
        is_fwd(true),
        kmer(length), current(length) {}

    Sequence(const ModelType &model_, i32 ref_id_, IntervalIndex<i64> coords_, bool is_fwd_) : 
        model(model_), 
        ref_id(ref_id_), 
        coords(coords_),
        is_fwd(is_fwd_), 
        kmer(coords.length), 
        current(coords.length) {}

    Sequence(const ModelType &model_, const std::string &seq) :
            Sequence(model_, seq.size()-K+1) {

        kmer[0] = model.str_to_kmer(seq);

        for (size_t i = 0; i < size()-1; i++) {
            auto b = BASE_BYTES[seq[i+ModelType::KMER_LEN]];
            kmer[i+1] = model.kmer_neighbor(kmer[i], b);
        }
        init_current();
    }

    //TODO input pacseq and interval index, set from each segment
    //then RefIndex can just feed right in
    //need to figure out mrefs, k-mer trim
    
    //eventually need to write new FastaIndex based on FAI
    //then also pass file pointer and interval index, read chunks

    Sequence(const ModelType &model_, u8 *seq, size_t start, size_t end) :
            Sequence(model_, end-start-ModelType::KMER_LEN+1) {
        kmer = model.pacseq_to_kmers(seq, start, end);
    }

    typename ModelType::KmerTypePy kmer_to_str(size_t i) const {
        return model.kmer_to_arr(i);
    }

    KmerType get_kmer(i64 r) const {
        auto i = coords.get_index(r);
        return kmer[i];
    }

    KmerType get_current(i64 r) const {
        auto i = coords.get_index(r);
        return current[i];
    }

    void init_current() {
        for (size_t i = 0; i < size(); i++) {
            current[i] = model.kmer_current(kmer[i]);
        }
    }

    size_t size() const {
        return kmer.size();
    }

    static void pybind(py::module &m, std::string suffix) {
        py::class_<Sequence> c(m, ("Sequence"+suffix).c_str());
        //auto c = super::template pybind<Sequence>(m, ("Sequence"+suffix).c_str(), false);

        c.def(py::init<const ModelType &, const std::string &>());
        c.def("__len__", &Sequence::size);
        c.def_readonly("coords", &Sequence::coords);
        c.def_readonly("ref_id", &Sequence::ref_id);
        c.def_readonly("kmer", &Sequence::kmer);
        c.def_readonly("current", &Sequence::current);
        c.def_readonly("is_fwd", &Sequence::is_fwd);
        c.def("kmer_to_str", py::vectorize(&Sequence::kmer_to_str));
        c.def("get_kmer", py::vectorize(&Sequence::get_kmer));
        c.def("get_current", py::vectorize(&Sequence::get_current));
        c.attr("K") = pybind11::cast(Sequence::K);
    }
};

#endif
