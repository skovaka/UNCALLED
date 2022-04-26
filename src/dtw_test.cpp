#include <iostream>
#include <math.h>
#include <unordered_map>

#include "event_profiler.hpp"
#include "normalizer.hpp"
#include "fast5_reader.hpp"
#include "ref_index.hpp"
#include "dtw.hpp"

//const std::string CONF_DIR(std::getenv("UNCALLED_CONF")),
//                  DEF_MODEL = CONF_DIR + "/r94_5mers.txt",
//                  DEF_CONF = CONF_DIR + "/defaults.toml";
//
//bool load_conf(int argc, char** argv, Config &config);


typedef struct {
    std::string rd_name, rf_name;
    u32 rd_st, rd_en;
    u64 rf_st, rf_en;
    bool fwd;
} Query;

std::unordered_map<std::string, Query> 
    load_queries(const std::string &fname, 
                 Fast5Iter &fast5s) {

    std::unordered_map<std::string, Query> queries;

    if (fname.empty()) {
        std::cerr << "Must specify query file\n";
        return queries;
    }
    
    std::ifstream infile(fname);

    if (!infile.is_open()) {
        std::cerr << "Error: failed to open query file\n";
        return queries;
    }

    Query q;
    char strand;
    while (!infile.eof()) {
        infile >> q.rd_name >> q.rd_st >> q.rd_en
               >> q.rf_name >> q.rf_st >> q.rf_en
               >> strand;
        q.fwd = strand == '+';

        queries[q.rd_name] = q;
        fast5s.add_read(q.rd_name);
    }

    return queries;
}

int main(int argc, char** argv) {
    Timer t;

    std::string index_prefix  = std::string(argv[1]), 
                fast5_fname   = std::string(argv[2]),
                query_fname = std::string(argv[3]);

    std::string out_prefix = "";
    if (argc > 4) {
        out_prefix = std::string(argv[4]);
    }

    Fast5Dict fast5s_test;
    fast5s_test.load_index(fast5_fname);
    auto rd = fast5s_test[query_fname];

    //if (false) {

    auto model = PoreModel<KmerLen::k5>("r94_dna");

    EventDetector evdt;
    EventProfiler evpr;

    RefIndex<KmerLen::k5> idx(index_prefix);
    idx.load_pacseq();

    Fast5Iter fast5s;
    fast5s.add_fast5(fast5_fname);

    bool create_events = true;

    auto queries = load_queries(query_fname, fast5s);

    while (!fast5s.empty()) {
        //Get next read and corrasponding query
        auto read = fast5s.next_read();
        std::cerr << read.get_id() << "\n";
        //std::cerr << "aligning " << read.get_id() << "\n";
        //std::cerr.flush();

        Query q = queries[read.get_id()];

        std::vector<u16> kmers = idx.get_kmers(q.rf_name, q.rf_st, q.rf_en);
        if (!q.fwd) kmers = kmers_revcomp<KmerLen::k5>(kmers);

        float read_mean = 0;
        for (u16 k : kmers) {
            read_mean += model.kmer_current(k);
        }
        read_mean /= kmers.size();
        float read_stdv = 0;
        for (u16 k : kmers) {
            read_stdv += pow(model.kmer_current(k) - read_mean, 2);
        }
        read_stdv = sqrt(read_stdv / kmers.size());
        Normalizer norm(read_mean, read_stdv);

        //Normalizer norm(model.model_mean(), model.model_stdv());

        //Get raw signal
        auto &full_raw = read.get_signal();
        std::vector<float> signal;
        if (q.rd_st != 0 || q.rd_en != 0) {
            u32 en = q.rd_en == 0 ? read.size() : q.rd_en;
            signal.reserve(en - q.rd_st);

            for (u32 i = q.rd_st; i < en; i++) {
                signal.push_back(full_raw[i]);
            }
        } else {
            signal = full_raw;
        }

        //Create events if needed
        if (create_events) {
            auto events = evdt.get_events(signal);
            auto mask = evpr.get_full_mask(events);
            signal.clear();

            for (u32 i = 0; i < events.size(); i++) {
                if (mask[i]) signal.push_back(events[i].mean);
            }
            
            norm.set_signal(signal);

        } else {
            norm.set_signal(signal);
        }

        //Normalize
        signal.clear();
        signal.reserve(norm.unread_size());
        while (!norm.empty()) signal.push_back(norm.pop());

        //Takes up too much space :(
        if (signal.size() > 50000) {
            //std::cerr << "Skipping " << read.get_id() << "\n";
            continue;
        }

        auto band_width = static_cast<i32>(signal.size() * 0.05);

        auto prms = DTW_PRMS_DEF;
        prms.band_width = 500;

        //StaticBDTW dtw(prms, signal, kmers, model);
        DTWd dtw(signal, kmers, model, prms);

        if (!out_prefix.empty()) {
            std::string path_fname = out_prefix+read.get_id()+".txt";
            std::ofstream out(path_fname);
            for (auto &t : dtw.get_path()) {
                out << t.ref << "\t" << t.qry << "\n";
            }
            out.close();
        }

        std::cout << read.get_id() << "\t"
                  << dtw.mean_score() << "\t"
                  << dtw.path_.size() << "\t"
                  << (t.lap()/1000) << "\n";
        std::cout.flush();


    }


    //}//end if false

    return 0;
}
