#ifndef EVENTALIGN_INCL
#define EVENTALIGN_INCL

#include "config.hpp"
#include "pore_model.hpp"
#include "intervals.hpp"

//template <KmerLen K, typename KmerType=typename std::conditional<(K < 8), u16, u32>::type>
//struct Eventalign : public DataFrame<i64, KmerType, i32, i32, float> {
//    using Super = DataFrame<i64, KmerType, i32, i32, float>;
//
//    template <size_t I>
//    using ColType = Super::typename ColType<I>;// = typename Super::ColType<I>;
//
//    static constexpr typename Super::NameArray names = {"pac", "model_kmer", "event_index", "standardized_level"}; 
//    ColType<0> &pac = std::get<0>(Super::data_); 
//    ColType<1> &model_kmer = std::get<1>(Super::data_);                      
//    ColType<2> &event_index = std::get<2>(Super::data_);                      
//    ColType<3> &standardized_level = std::get<3>(Super::data_);
//    using Super::DataFrame;                              
//};                
//constexpr Eventalign::NameArray Eventalign::names;

//void write_eventalign(Eventalign<K> df, ProcessedRead read, i32 read_idx, bool fwd) {
template <typename ModelType>
std::string write_eventalign_new(Alignment<ModelType> &aln, bool write_name, bool signal_index, py::array_t<float> signal_np) {

    auto signal = PyArray<float>(signal_np);
    std::stringstream ss;

    auto coord = aln.seq.coord;
    auto &model = aln.seq.model;
    auto sample_rate = model.PRMS.sample_rate;

    for (size_t i = 0; i < aln.dtw.size(); i++) {
        auto &samps = aln.dtw.samples.coords[i];
        if (!samps.is_valid()) {
            continue;
        }
        //if (aln.dtw.current[i] == ValArray<float>::NA) continue;
        //auto &evt = read.events[i];
        auto kmer = aln.seq.kmer[i], model_kmer = kmer;
        if (model.PRMS.reverse) model_kmer = model.kmer_rev(kmer);
        auto ref_kmer = model_kmer;
        if (!coord.fwd()) ref_kmer = model.kmer_revcomp(ref_kmer);

        auto ref = aln.seq.mpos[i];
        if (ref < 0) ref = -ref-1;


        ss << coord.name << "\t"
           << ref - model.PRMS.shift << "\t"
           << model.kmer_to_str(ref_kmer) << "\t";

        if (write_name) {
           ss << aln.read_id;
        } else {
           ss << aln.id;
        }
        ss << "\tt" << "\t"
           << i << "\t"
           << model.current.norm_to_pa(aln.dtw.current[i]) << "\t"
           << (model.current.norm_to_pa_sd(aln.dtw.current_sd[i])) << "\t"
           << (samps.length() / sample_rate) << "\t"
           << model.kmer_to_str(model_kmer) << "\t"
           << model.current.norm_to_pa(model.current.mean[kmer]) << "\t"
           << model.current.norm_to_pa_sd(model.current.stdv[kmer]) << "\t"
           << aln.dtw.current[i];

        if (signal_index) {
           ss << "\t" << samps.start << "\t" << samps.end; //<< "\n";
        }

        if (signal.size() > 0) {
            ss << "\t" << signal[samps.start];
            for (size_t j = samps.start+1; j < samps.end; j++) {
                ss << "," << signal[j];
            }
        }

        ss << "\n";
    }
    return ss.str();
}

template <typename ModelType>
std::string write_eventalign(
        Config &conf,
        ModelType &model,
        std::string read_id,
        bool fwd,
        ProcessedRead read, 
        const std::string &ref_name, py::array_t<i64> ref_np,
        bool signal_index,
        py::array_t<typename ModelType::kmer_t> kmer_np,
        py::array_t<i32> event_index_np,
        py::array_t<float> std_level_np,
        py::array_t<float> signal_np) {

    auto ref = PyArray<i64>(ref_np);
    auto kmers = PyArray<typename ModelType::kmer_t>(kmer_np);
    auto event_index = PyArray<i32>(event_index_np);
    auto std_level = PyArray<float>(std_level_np);
    auto signal = PyArray<float>(signal_np);

    float sample_rate = conf.pore_model.sample_rate;

    std::stringstream ss;

    for (size_t i = 0; i < ref.size(); i++) {
        auto &evt = read.events[i];
        auto kmer = kmers[i], model_kmer = kmer;
        if (conf.pore_model.reverse) model_kmer = model.kmer_rev(kmer);
        auto ref_kmer = model_kmer;
        if (!fwd) ref_kmer = model.kmer_comp(ref_kmer);

        ss << ref_name << "\t"
           << ref[i] << "\t"
           << model.kmer_to_str(ref_kmer) << "\t"
           << read_id << "\t"
           << "t" << "\t"
           << event_index[i] << "\t"
           << model.current.norm_to_pa(evt.mean) << "\t"
           << (model.current.norm_to_pa_sd(evt.stdv)) << "\t"
           << (evt.length / sample_rate) << "\t"
           << model.kmer_to_str(model_kmer) << "\t"
           << model.current.norm_to_pa(model.current.mean[kmer]) << "\t"
           << model.current.norm_to_pa_sd(model.current.stdv[kmer]) << "\t"
           << std_level[i];

        if (signal_index) {
           ss << "\t" << evt.start << "\t"
              << evt.start + evt.length; //<< "\n";
        }

        if (signal.size() > 0) {
            ss << "\t" << signal[evt.start];
            for (size_t j = 1; j < evt.length; j++) {
                ss << "," << signal[evt.start+j];
            }
        }

        ss << "\n";
    }
    return ss.str();
}

#endif
