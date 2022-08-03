#ifndef EVENTALIGN_INCL
#define EVENTALIGN_INCL

#include "config.hpp"
#include "pore_model.hpp"
#include "dataframe.hpp"

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
std::string write_eventalign(
        Config &conf,
        ModelType &model,
        i32 read_idx, bool fwd,
        ProcessedRead read, 
        const std::string &ref_name, py::array_t<i64> ref_np,
        py::array_t<typename ModelType::kmer_t> kmer_np,
        py::array_t<i32> event_index_np,
        py::array_t<float> std_level_np,
        std::string read_id,
        py::array_t<float> signal_np) {

    auto ref = PyArray<i64>(ref_np);
    auto kmers = PyArray<typename ModelType::kmer_t>(kmer_np);
    auto event_index = PyArray<i32>(event_index_np);
    auto std_level = PyArray<float>(std_level_np);
    auto signal = PyArray<float>(signal_np);

    float sample_rate = conf.read_buffer.sample_rate;

    std::stringstream ss;

    for (size_t i = 0; i < ref.size(); i++) {
        auto &evt = read.events[i];
        auto kmer = kmers[i], model_kmer = kmer;
        if (!conf.read_buffer.seq_fwd) model_kmer = model.kmer_rev(kmer);
        auto ref_kmer = model_kmer;
        if (!fwd) ref_kmer = model.kmer_comp(ref_kmer);

        ss << ref_name << "\t"
           << ref[i] << "\t"
           << model.kmer_to_str(ref_kmer) << "\t";

        if (read_id.empty()) {
            ss << read_idx; 
        } else {
            ss << read_id;
        }

        ss  << "\t" << "t" << "\t"
            << event_index[i] << "\t"
            << evt.mean << "\t"
            << evt.stdv << "\t"
            << (evt.length / sample_rate) << "\t"
            << ModelType::kmer_to_str(model_kmer) << "\t"
            << model.kmer_means_[kmer] << "\t"
            << model.kmer_stdvs_[kmer] << "\t"
            << std_level[i];
            //<< evt.start << "\t"
            //<< evt.start + evt.length; //<< "\n";

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

//"contig" : track.coords.ref_name,
//"position" : events.index-2,
//"reference_kmer" : ref_kmers,
//"read_index" : self.prev_aln_id,
//"strand" : "t",
//"event_index" : pd.RangeIndex(0,len(events))[::-1]+1,
//"event_level_mean" : events["current"],
//"event_stdv" : stdvs,
//"event_length" : events["length"] / track.conf.read_buffer.sample_rate,
//"model_kmer" : model_kmers,
//"model_mean" : model.means[kmers],
//"model_stdv" : model.stdvs[kmers],
//"standardized_level" : std_level,
//"start_idx" : events["start"],
//"end_idx" : events["start"] + events["length"],

#endif
