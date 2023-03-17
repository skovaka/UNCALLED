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
#include <numeric>
#include "util.hpp"
#include "intervals.hpp"

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
    int k, shift;
    float pa_mean, pa_stdv, norm_max;
    bool reverse, complement;

    template<typename T>
    T norm_to_pa(T norm) const {
        return (norm * pa_stdv) + pa_mean;
    }

    template<typename T>
    T pa_to_norm(T val) const {
        return (val - pa_mean) / pa_stdv;
    }

    template<typename T>
    T pa_sd_to_norm(T val) const {
        return val / pa_stdv;
    }

    template<typename T>
    T norm_to_pa_sd(T val) const {
        return val * pa_stdv;
    }
};

extern const PoreModelParams PORE_MODEL_PRMS_DEF;

struct ModelDF {
    using inorm_t = i16;
    static constexpr inorm_t INORM_MAX = std::numeric_limits<inorm_t>::max(),
                             INORM_NA = std::numeric_limits<inorm_t>::min();
    
    std::reference_wrapper<PoreModelParams> PRMS;
    float INORM_SCALE;
    ValArray<float> mean, stdv;

    static float init_inorm(const PoreModelParams &prms) {
        return INORM_MAX / prms.norm_max;
    }

    //ModelDF(const &ModelDF df) = default; 

    ModelDF(PoreModelParams &params, size_t len=0, float fill=0) : 
        PRMS(params),
        INORM_SCALE(init_inorm(params)),
        mean(fill, len), stdv(fill, len) {
        //mean_shift(0), mean_scale(1) {
    }

    template <typename Container>
    ModelDF(PoreModelParams &params, std::vector<size_t> &order, Container &means) :
            PRMS(params),
            INORM_SCALE(init_inorm(params)) {
        set_means(order, means);
    }

    //ModelDF(const PoreModelParams &params, PyArray<float> means, PyArray<float> stdvs, bool normalize) {
    template <typename Container>
    ModelDF(PoreModelParams &params, std::vector<size_t> &order, Container &means, Container &stdvs) :
            PRMS(params),
            INORM_SCALE(init_inorm(params)) {
        //PyArray<float> means_py(means), stdvs_py(stdvs);
        set_means(order, means);
        set_stdvs(order, stdvs);
    }

    //std::pair<float, float> get_norm_params() const {
    //    auto shift = mean.mean(),
    //         scale = mean.stdv(shift);
    //    return {scale, shift};
    //}

    void normalize() {
        mean = pa_to_norm(mean);
        stdv = pa_sd_to_norm(stdv);
    }
    
    template <typename Container>
    void set_means(std::vector<size_t> &order, Container &vals) {
        if (vals.size() == 0) return;
        set_vals(order, mean, vals);
        //if (normalize && mean.size() > 0) {
        //    mean_shift = mean.mean(); 
        //    mean_scale = mean.stdv(mean_shift);
        //    mean = pa_to_norm(mean);
        //} else {
        //    mean_shift = 0;
        //    mean_scale = 1;
        //}
    }

    template <typename Container>
    void set_stdvs(std::vector<size_t> &order, Container &vals) {
        if (vals.size() == 0) return;
        set_vals(order, stdv, vals);
       // if (normalize && stdv.size() > 0) {
       //     stdv /= mean_scale;
       // }
    }

    template <typename Container>
    void set_vals(std::vector<size_t> &order, ValArray<float> &dest, const Container &vals) {
        if (dest.size() > 0 && dest.size() != vals.size()) {
            throw std::runtime_error("Invalid PoreModel length");
        }
        if (order.size() == 0) {
            dest = ValArray<float>(vals.begin(), vals.size());
        } else if (order.size() != vals.size()) {
            throw std::runtime_error("Order and values mismatch");
        } else {
            dest =  ValArray<float>(vals.size());
            size_t i = 0;
            for (auto j : order) {
                dest[i++] = vals[j];
            }
        }
    }

    float get_mean_pa(size_t i) {
        return norm_to_pa(mean[i]);
    }

    float get_stdv_pa(size_t i) {
        return norm_to_pa(stdv[i]);
    }

    inorm_t get_mean_inorm(size_t i) {
        return norm_to_inorm(mean[i]);
    }

    inorm_t get_stdv_inorm(size_t i) {
        return norm_to_inorm(stdv[i]);
    }

    template<typename T>
    T norm_to_pa(T val) const {
        return PRMS.get().norm_to_pa(val);
        //return (norm * PRMS.get().pa_stdv) + PRMS.get().pa_mean;
    }

    template<typename T>
    T pa_to_norm(T val) const {
        return PRMS.get().pa_to_norm(val);
        //(val - PRMS.get().pa_mean) / PRMS.get().pa_stdv;
    }

    template<typename T>
    T pa_sd_to_norm(T val) const {
        return PRMS.get().pa_sd_to_norm(val);
        //return val / PRMS.get().pa_stdv;
    }

    template<typename T>
    T norm_to_pa_sd(T val) const {
        return PRMS.get().norm_to_pa_sd(val);
        //return val * PRMS.get().pa_stdv;
    }

    inorm_t norm_to_inorm(float norm) const {
        return static_cast<inorm_t>(round(norm * INORM_SCALE));
    }

    float inorm_to_norm(inorm_t norm) const {
        return norm / INORM_SCALE; //(norm_max / INORM_MAX);
    }

    size_t size() const {
        return mean.size();
    }

    static void pybind(py::module_ &m) {
        py::class_<ModelDF> c(m, "ModelDF");
        c.def(py::init<PoreModelParams &, size_t, float>());
        c.def(py::init<PoreModelParams &, size_t, float>());
        c.def(py::init<PoreModelParams &, size_t>());

        //c.def(py::init<const PoreModelParams &, py::array_t<float>, py::array_t<float>, bool>());
        //c.def(py::init<const PoreModelParams &, py::array_t<float>, bool>());


        //c.def(py::init<>());
        c.def_readonly("mean", &ModelDF::mean);
        c.def_readonly("stdv", &ModelDF::stdv);
        //c.def_readonly("mean_scale", &ModelDF::mean_scale);
        //c.def_readonly("mean_shift", &ModelDF::mean_shift);
        //c.def_readonly("stdv_scale", &ModelDF::stdv_scale);
        //c.def("set_means", static_cast< void (*) (ModelDF::set_means >);
        //        py::vectorize(static_cast< KmerType (*) (const KmerTypePy &, u32)>(&Class::str_to_kmer)), 
        //c.def("set_stdvs", &ModelDF::set_stdvs);
        c.def_readonly("INORM_SCALE", &ModelDF::INORM_SCALE);
        c.def("norm_to_inorm", py::vectorize(&ModelDF::norm_to_inorm));
        c.def("get_mean_inorm", py::vectorize(&ModelDF::get_mean_inorm));
        c.def("get_mean_pa", py::vectorize(&ModelDF::get_mean_pa));
        c.def("get_stdv_pa", py::vectorize(&ModelDF::get_stdv_pa));
        c.def("__len__", &ModelDF::size);
    }
};

//template<KmerLen KMER_LEN, typename KmerType=typename std::conditional<(KMER_LEN < 8), u16, u32>::type>
template<typename KmerType>
class PoreModel {

    public:
    const KmerLen KMER_LEN;
    const KmerType KMER_MASK, KMER_COUNT;

    using kmer_t = KmerType;

    PoreModelParams PRMS;

    ModelDF current;

    //current.mean, current.stdv, 
    ValArray<KmerType> kmer;
    std::vector<float> current_var2x, lognorm_denoms;
    float model_mean_, model_stdv_;
    bool loaded_, compl_;

    PoreModel(PoreModelParams p) : 
        KMER_LEN(p.k),
        KMER_MASK((1 << (2*KMER_LEN)) - 1),
        KMER_COUNT(static_cast<KmerType>(pow(BASE_COUNT, KMER_LEN))),
        PRMS(p) ,
        current(PRMS, KMER_COUNT)
            {
        loaded_ = false;

        model_mean_ = 0;
        model_stdv_ = 1;

        //current.mean.resize(KMER_COUNT);
        //current.stdv.resize(KMER_COUNT);
        kmer.resize(KMER_COUNT);
        current_var2x.resize(KMER_COUNT);
        lognorm_denoms.resize(KMER_COUNT);

        //if (p.name.empty()) return;

        //load_presets();
        //auto preset = PORE_MODEL_PRESETS.find(PRMS.name);

        //if (preset != PORE_MODEL_PRESETS.end()) {
        //    init_vals(preset->second.vals);
        //} else {
        //init_tsv(PRMS.name);
        //}
    }

    template <typename Fn>
    std::vector<size_t> kmer_order(size_t n, Fn kmer_fn) {
        std::vector<size_t> order(n);
        std::iota(order.begin(), order.end(), 0);
        
        std::stable_sort(order.begin(), order.end(), 
            [&kmer_fn](size_t a, size_t b) {return kmer_fn(a) < kmer_fn(b);});
        return order;
    }

    PoreModel(const PoreModel &model) = default; 

    PoreModel(PoreModelParams p, py::array_t<float> current_mean_py, py::array_t<float> current_stdv_py, bool normalize) : PoreModel(p) {
        PyArray<float> current_mean(current_mean_py), current_stdv(current_stdv_py);

        std::vector<size_t> order; //(current_mean.size());
        //std::iota(order.begin(), order.end(), 0);
        auto n = current_mean.size();
        auto get_kmer = [&](size_t i) {return i;};
        if (PRMS.reverse && PRMS.complement) {
            order = kmer_order(n, [&](size_t i) {return kmer_revcomp(i);});
        } else if (PRMS.reverse) {
            order = kmer_order(n, [&](size_t i) {return kmer_rev(i);});
        } else if (PRMS.complement) {
            order = kmer_order(n, [&](size_t i) {return kmer_comp(i);});
        } else {
            order = kmer_order(n, [&](size_t i) {return i;});
        }

        current = ModelDF(PRMS, order, current_mean, current_stdv);

        if (normalize) {
            PRMS.pa_mean = current.mean.mean();
            PRMS.pa_stdv = current.mean.stdv(PRMS.pa_mean);
            current.normalize();
        }
    }

    PoreModel() : PoreModel(PORE_MODEL_PRMS_DEF) {
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

    /* Serializes pore model into the following format:
    u8  kmer_bytes: number of bytes per k-mer (kmer dtype)
    u32 kmer_count: number of kmers in model
    bool current: true if current is stored
    bool current_sd: true if current_sd is stored
    bool dwell: true if dwell is stored
    bool stdv: true if stdv is stored for current, current_sd, and dwell
    */
    //ValArray<u8> serialize() const {
    //    size_t size = sizeof(i32) * 
    //}
        
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

    float norm_to_pa(float val) const {
        return PRMS.norm_to_pa(val);
    }

    float pa_to_norm(float val) const {
        return PRMS.pa_to_norm(val);
    }

    float pa_sd_to_norm(float val) const {
        return PRMS.pa_sd_to_norm(val);
    }

    float norm_to_pa_sd(float val) const {
        return PRMS.norm_to_pa_sd(val);
    }

    float norm_pdf(float samp, KmerType kmer) const {
        return (-pow(samp - current.mean[kmer], 2) / current_var2x[kmer]) - lognorm_denoms[kmer];
    }

    float abs_diff(float samp, KmerType kmer) const {
        return std::abs(samp - current.mean[kmer]);
    }

    float z_score(float samp, KmerType kmer) const {
        return std::abs(samp - current.mean[kmer]) / current.stdv[kmer];
    }

    float model_mean() const {
        return model_mean_;
    }

    float model_stdv() const {
        return model_stdv_;
    }

    float kmer_current(KmerType kmer) const {
        return current.mean[kmer];
    }

    float kmer_stdv(KmerType kmer) const {
        return current.stdv[kmer];
    }

    bool is_loaded() const {
        return loaded_;
    }

    u32 get_kmer_count() const {
        return KMER_COUNT;
    }

    void init_stdv() {
        model_stdv_ = 0;

        for (KmerType kmer = 0; kmer < KMER_COUNT; kmer++) {
            model_stdv_ += pow(current.mean[kmer] - model_mean_, 2);
        }

        model_stdv_ = sqrt(model_stdv_ / KMER_COUNT);
    }

    void init_kmer(KmerType k, float mean, float stdv) {
        if (PRMS.reverse)  k = kmer_rev(k);
        if (PRMS.complement) k = kmer_comp(k);

        kmer[k] = k;
        current.mean[k] = mean;
        current.stdv[k] = stdv;
        current_var2x[k] = 2 * stdv * stdv;
        lognorm_denoms[k] = log(sqrt(M_PI * current_var2x[k]));
    }

    char base_to_char(u8 base) const {
        return BASE_CHARS[base];
    }

    u8 base_comp(u8 base) const {
        return BASE_COMP_B[base];
    }

    KmerType str_to_kmer(const std::string &kmer, u32 offs=0) const {
        KmerType index = BASE_BYTES[(u8) kmer[offs]];
        for (u8 i = 1; i < (u8) KMER_LEN; i++) {
            index = (index << 2) | BASE_BYTES[(u8) kmer[offs+i]];
        }
        return index;
    }

    KmerType kmer_comp(KmerType kmer) const {
        return kmer ^ KMER_MASK;
    }

    private:
    u16 _kmer_rev(u16 kmer) const {
        kmer = ( (kmer >> 2 & 0x3333) | (kmer & 0x3333) << 2 );
        kmer = ( (kmer >> 4 & 0x0F0F) | (kmer & 0x0F0F) << 4 );
        kmer = ( (kmer >> 8 & 0x00FF) | (kmer & 0x00FF) << 8 );
        return kmer >> (2 * (8 - KMER_LEN));
    }

    u32 _kmer_rev(u32 kmer) const {
        kmer = ( (kmer >> 2  & 0x33333333) | (kmer & 0x33333333) << 2 );
        kmer = ( (kmer >> 4  & 0x0F0F0F0F) | (kmer & 0x0F0F0F0F) << 4 );
        kmer = ( (kmer >> 8  & 0x00FF00FF) | (kmer & 0x00FF00FF) << 8 );
        kmer = ( (kmer >> 16 & 0x0000FFFF) | (kmer & 0x0000FFFF) << 16 );
        return kmer >> (2 * (16 - KMER_LEN));
    }

    public:
    KmerType kmer_rev(KmerType kmer) const {
        return _kmer_rev(kmer);
    }

    KmerType kmer_revcomp(KmerType kmer) const {
        return kmer_rev(~kmer);
    }

    std::vector<KmerType> kmers_revcomp(const std::vector<KmerType> &kmers) const {
        std::vector<KmerType> rev;
        rev.reserve(kmers.size());
        for (auto k = kmers.rbegin(); k != kmers.rend(); k++) {
            rev.push_back(kmer_revcomp(*k));
        }
        return rev;
    }

    u8 kmer_head(KmerType kmer) const {
        return (u8) ((kmer >> (2*( (u8) KMER_LEN ) - 2)) & 0x3);
    }

    KmerType kmer_neighbor(KmerType kmer, u8 i) const {
        return ((kmer << 2) & KMER_MASK) | i; 
    }

    KmerType set_kmer_base(KmerType kmer, u8 i, u8 b) const {
        auto shift = 2*i;
        auto mask = ~(0x3 << shift);
        return (kmer & mask) | (b << shift); 
    }

    u8 kmer_base(KmerType kmer, u8 i) const {
        return (u8) ((kmer >> (2 * ((KmerType)KMER_LEN-i-1))) & 0x3);
    }

    u8 kmer_base_count(KmerType kmer, u8 base) const {
        u8 count = 0;
        for (auto i = 0; i < KMER_LEN; i++) {
            count += kmer_base(kmer, i) == base;
        }
        return count;
    }

    std::string kmer_to_str(KmerType kmer) const {
        std::string s(KMER_LEN, 'N');
        for (u8 i = 0; i < KMER_LEN; i++) {
            s[i] = BASE_CHARS[kmer_base(kmer, i)];
        }
        return s;
    }

    ValArray<KmerType> pacseq_to_kmers(u8 *seq, u64 st, u64 en) const {
        ValArray<KmerType> ret(en - st - KMER_LEN + 1);

        u64 pst = st >> 2,
            pen = ((en) >> 2)+1;

        u64 i = 0;
        KmerType kmer = 0;
        u8 bst = (st&3), ben;

        //TODO could optimize by splitting into init and extend loops
        for (u64 j = pst; j < pen; j++) {
            ben = j == pen-1 ? (en&3) : 4;
            for (u8 k = bst; k < ben; k++) {
                kmer = kmer_neighbor(kmer, (seq[j] >> ((KMER_LEN^3) << 1) ) & 3);
                if (i >= KMER_LEN-1) {
                    ret[i] = kmer;
                }
                i++;
            }
            bst = 0;
        }

        return ret;
    }

    #ifdef PYBIND

    #define PY_MODEL_DEF(P, D) c.def(#P, &PoreModel::P, D);
    #define PY_MODEL_DEF_STATIC(P, D) c.def_readwrite_static(#P, &PoreModel::P, D);
    #define PY_MODEL_DEF_CONSTANT(P, D) c.def_readonly_static(#P, &PoreModel::P, D);
    #define PY_MODEL_DEFVEC(P, D) c.def(#P, py::vectorize(&PoreModel::P), D);
    #define PY_MODEL_PROP(P, D) c.def_property_readonly(#P, &PoreModel::P, D);
    #define PY_MODEL_PARAM(P, D) p.def_readwrite(#P, &PoreModelParams::P, D);

    public:

	//using KmerTypePy = std::array<char, KMER_LEN+1>;

	//KmerTypePy kmer_to_arr(KmerType kmer) {
	//	KmerTypePy ret{0};
	//	for (size_t i = 0; i < KMER_LEN; i++) {
	//		ret[i] = BASE_CHARS[kmer_base(kmer, i)];
	//	}
	//	return ret;
    //}

    //KmerType str_to_kmer(const KmerTypePy &kmer, u32 offs) {
	//	auto s = std::string(kmer.data());
	//	return str_to_kmer(s, offs);
	//}



    static void pybind(py::module_ &m, const std::string &suffix) {
        using Class = PoreModel<KmerType>;

        py::class_<Class> c(m, ("PoreModel" + suffix).c_str());

        c.def(pybind11::init<const Class &>());
        //c.def(pybind11::init<KmerLen>());
        c.def(pybind11::init<PoreModelParams>());
        c.def(pybind11::init<PoreModelParams, py::array_t<float>, py::array_t<float>, bool>());
        //c.def(pybind11::init<const std::string &, bool, bool, KmerLen>(), 
        //      py::arg("name"), py::arg("reverse")=PORE_MODEL_PRMS_DEF.reverse, py::arg("complement")=PORE_MODEL_PRMS_DEF.complement, py::arg("k")=PORE_MODEL_PRMS_DEF.k);
        //c.def(pybind11::init<const std::vector<float> &, bool, bool, KmerLen>(), 
        //      py::arg("vals"), py::arg("reverse")=PORE_MODEL_PRMS_DEF.reverse, py::arg("complement")=PORE_MODEL_PRMS_DEF.complement, py::arg("k")=PORE_MODEL_PRMS_DEF.k);
        
        //c.attr("KmerArray") = m.attr("PyArray;

        c.def_property_readonly("kmer_count", &Class::get_kmer_count, "The number of k-mers in the model");
        c.def_readwrite("PRMS", &Class::PRMS, "Class parameters");
        PY_MODEL_PROP(model_mean, "The mean of all model k-mer currents");
        PY_MODEL_PROP(model_stdv, "The standard deviation of all model k-mer currents");

        //c.def_static("is_preset", &Class::is_preset, "List of model preset names");
        //c.def_static("get_preset_names", &Class::get_preset_names, "List of model preset names");

        //c.def_property_readonly("means", 
        //    [](Class &r) -> pybind11::array_t<float> {
        //        return pybind11::array_t<float>(r.current.mean.size(), r.current.mean.data());
        //}, "The expected mean current of each k-mer");

        //c.def_property_readonly("stdvs",
        //    [](Class &r) -> pybind11::array_t<float> {
        //        return pybind11::array_t<float>(r.current.stdv.size(), r.current.stdv.data());
        //}, "The expected standard devaition of each k-mer");

        c.def_readonly("current", &PoreModel::current);

        PY_MODEL_DEFVEC(norm_pdf, "Returns the log probability that the current matches the k-mer based on the normal distibution probability density function");

        PY_MODEL_DEFVEC(abs_diff, "Returns the absolute difference between the observed and model current");

        c.def("__getitem__", py::vectorize(&Class::kmer_current), "Alias for get_current()");

        c.def("__len__", &Class::get_kmer_count, "Alias for kmer_count");

        //c.def_readonly_static("KMER_LEN", &Class::KMER_LEN);
        //c.def_readonly_static("KMER_COUNT", &Class::KMER_COUNT);
        c.def_readonly("K", &Class::KMER_LEN);
        c.def_readonly("KMER_COUNT", &Class::KMER_COUNT);

        c.def("str_to_kmer", &PoreModel::str_to_kmer);
        //        py::vectorize(static_cast< KmerType (*) (const KmerTypePy &, u32)>(&Class::str_to_kmer)), 
        //        py::arg("kmer"), py::arg("offs")=0);

        c.def("kmer_to_str",  &Class::kmer_to_str, "Convert binary k-mers to strings");
        //c.def("kmer_to_arr",  py::vectorize(&Class::kmer_to_arr));

        c.def("base_to_char", py::vectorize(&Class::base_to_char), "Convert base index to base character"); 
        c.def("base_comp", py::vectorize(&Class::base_comp), "Returns the complement of a base index");
        
        c.def("kmer_rev",      py::vectorize(&Class::kmer_rev), "Reverses binary k-mers (not complement)");
        c.def("kmer_comp",     py::vectorize(&Class::kmer_comp), "Complements binary k-mers (not reverse)");
        c.def("kmer_revcomp",  py::vectorize(&Class::kmer_revcomp), "Reverse complements binary k-mers");
        c.def("kmer_head",     py::vectorize(&Class::kmer_head), "Returns the first base index in binary k-mers");
        c.def("kmer_base",     py::vectorize(&Class::kmer_base), "Returns the base at the specified index in a binary k-mer");
        c.def("set_kmer_base",     py::vectorize(&Class::set_kmer_base), "Sets a binary k-mer base");
        c.def("kmer_base_count",     py::vectorize(&Class::kmer_base_count), "Returns the number of occurances of the specified base index in the binary k-mer");
        c.def("kmer_neighbor", py::vectorize(&Class::kmer_neighbor), "Returns the binary k-mer shifted left with the specified base appended");

        c.def("norm_to_pa",     py::vectorize(&Class::norm_to_pa), "Convert normalized current levels to picoamps");
        c.def("norm_to_pa_sd",     py::vectorize(&Class::norm_to_pa_sd), "Convert normalized current standard deviations to picoamps");

        c.def("pa_to_norm",     py::vectorize(&Class::pa_to_norm), "Convert pA current levels to 0-mean 1-stdv normalized");
        c.def("pa_sd_to_norm",     py::vectorize(&Class::pa_sd_to_norm), "Convert picoamps standard deviations to 0-mean 1-stdv normalized");

    }

    #endif
};

template <typename ModelType>
struct Sequence {//: public DataFrame<typename ModelType::kmer_t, float, u8> {
    using KmerType = typename ModelType::kmer_t;
    //using super = DataFrame<KmerType, u8, float>;

    const ModelType &model;
    const KmerType KMER_LEN; //= ModelType::KMER_LEN;
    i32 ref_id;
    IntervalIndex<i64> coords;
    bool is_fwd; //TODO infer from mpos coords

    //static constexpr typename super::NameArray names = {"ref", "start", "end"}; 
    //typename super::template ColType<0> &kmer = std::get<0>(super::data_);   
    //typename super::template ColType<1> &current = std::get<1>(super::data_);   
    //typename super::template ColType<2> &base = std::get<2>(super::data_);   
    //ValArray<u8> base;

    ValArray<KmerType> kmer; 
    ValArray<float> current;   

    Sequence(const ModelType &model_, size_t length) : 
        model(model_), 
        KMER_LEN(model_.KMER_LEN),
        ref_id(-1),
        coords({{0,static_cast<i64>(length)}}), 
        is_fwd(true),
        kmer(length), current(length) {}

    Sequence(const ModelType &model_, i32 ref_id_, IntervalIndex<i64> coords_, bool is_fwd_) : 
        model(model_), 
        KMER_LEN(model_.KMER_LEN),
        ref_id(ref_id_), 
        coords(coords_),
        is_fwd(is_fwd_), 
        kmer(coords.length), 
        current(coords.length) {}

    Sequence(const ModelType &model_, const std::string &seq) :
            Sequence(model_, seq.size()-KMER_LEN+1) {

        kmer[0] = model.str_to_kmer(seq);

        for (size_t i = 0; i < size()-1; i++) {
            auto b = BASE_BYTES[seq[i+KMER_LEN]];
            kmer[i+1] = model.kmer_neighbor(kmer[i], b);
        }
        init_current();
    }

    //TODO input pacseq and interval index, set from each segment
    //then RefIndex can just feed right in
    //need to figure out mposs, k-mer trim
    
    //eventually need to write new FastaIndex based on FAI
    //then also pass file pointer and interval index, read chunks

    Sequence(const ModelType &model_, u8 *seq, size_t start, size_t end) :
            Sequence(model_, end-start-KMER_LEN+1) {
        kmer = model.pacseq_to_kmers(seq, start, end);
    }

    //typename ModelType::KmerTypePy kmer_to_str(size_t i) const {
    //    return model.kmer_to_arr(i);
    //}

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

    const ModelType &get_model() {
        return model;
    }

    static void pybind(py::module &m, std::string suffix) {
        py::class_<Sequence> c(m, ("Sequence"+suffix).c_str());
        //auto c = super::template pybind<Sequence>(m, ("Sequence"+suffix).c_str(), false);

        c.def(py::init<const ModelType &, const std::string &>());
        c.def("__len__", &Sequence::size);
        c.def_property_readonly("model", &Sequence::get_model);
        c.def_readonly("K", &Sequence::KMER_LEN);
        c.def_readonly("coords", &Sequence::coords);
        c.def_readonly("ref_id", &Sequence::ref_id);
        c.def_readonly("kmer", &Sequence::kmer);
        c.def_readonly("current", &Sequence::current);
        c.def_readonly("is_fwd", &Sequence::is_fwd);
        //c.def("kmer_to_str", py::vectorize(&Sequence::kmer_to_str));
        c.def("get_kmer", py::vectorize(&Sequence::get_kmer));
        c.def("get_current", py::vectorize(&Sequence::get_current));
    }
};

#endif
