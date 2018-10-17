#include <iostream>
#include <string>
#include <list>
#include <cstdlib>
#include "arg_parse.hpp"
#include "bwa_fmi.hpp"
#include "timer.hpp"

void brute_align(BwaFMI &fmi, const std::vector<u8> &bases, u32 sample_dist) {
    for (u64 i = 0; i < bases.size(); i++) {
        if (rand() % sample_dist != 0) {
            continue;
        }

        Range r = fmi.get_full_range(bases[i]);
        u64 j = i+1;
        for (; j < bases.size() && r.length() > 1; j++) {
            std::cout << r.length() << "\t";
            r = fmi.get_neighbor(r, bases[j]);
        }
        std::cout << r.length() << std::endl;
    }
}

void dyn_align(BwaFMI &fmi, const std::vector<u8> &bases, u32 sample_dist) {
    std::list<Range> ranges;

    bool output = (rand() % sample_dist) == 0;

    //Align from first base until fm range reaches 1
    ranges.push_back(fmi.get_full_range(bases[0]));
    for (u64 i = 1; i < bases.size() && ranges.back().length() > 1; i++) {

        if (output) std::cout << ranges.back().length() << "\t";

        ranges.push_back( fmi.get_neighbor(ranges.back(), bases[i]) );
    }
    
    if (output) std::cout << ranges.back().length() << std::endl;

    //Find alignments starting from each base
    for (u64 i = 1; i < bases.size(); i++) {
        bool output = (rand() % sample_dist) == 0;

        ranges.pop_front();//No longer require previous first base
        
        Range l,r; //Will store ranges to fill in to left and right of prev
        r = fmi.get_full_range(bases[i]); //Get full range of this base
        l = r.split_range(ranges.front()); //Only keep the parts not covered

        u64 j = i + 1;//next base
        auto m = ranges.begin();//current middle range
        for (; j < bases.size() && (l.is_valid() || r.is_valid()) && m != ranges.end(); m++) {

            //Extend from left range
            if (l.is_valid()) {
                m->start_ = l.start_; //extend left edge of middle
                l = fmi.get_neighbor(l, bases[j]); //get next left range 
            }
    
            //Extend from right range
            if (r.is_valid()) {
                m->end_ = r.end_; //extend right edge of middle
                r = fmi.get_neighbor(r, bases[j]); //get next right range
            }

            //Output
            if (output) { 
                std::cout << m->length();
                if (m->length() > 1) std::cout << "\t";
            }

            j++;
        }

        //If m != end, right+left must have reached 0 before end
        //Means removing first base didn't add new alignments

        //Output previously found ranges again
        if (output) { 
            for (; m != ranges.end(); m++) {
                std::cout << m->length();
                if (m->length() > 1) std::cout << "\t";
            }
        }

        j = i + ranges.size();

        //Removing first base DID add new alignments
        while (ranges.back().length() > 1 && j < bases.size()) {
            r = fmi.get_neighbor(ranges.back(), bases[j]);
            if (!r.is_valid()) { 
                break;
            }
            ranges.push_back(r);
            if (output) { 
                std::cout << r.length();
                if (r.length() > 1) std::cout << "\t";
            }
            j++;
        }

        if (output) std::cout << std::endl;

    }
}

enum Opt {
    REF_FASTA = 'f',
    BWA_PREFIX = 'x',
    SAMPLE_DIST = 'r',
    RANDOM_SEED = 's',
    DYNAMIC_ALIGN = 'd'
};

int main(int argc, char** argv) {
    ArgParse args("Refererence self-aligner. Aligns a fasta reference to itself from random locations. For each starting location outputs the length of the FM index range at each step of the alignment until that length reaches 1 (a unique location is found).");
    args.add_string(Opt::REF_FASTA, "ref_fasta", "", "FASTA reference");
    args.add_string(Opt::BWA_PREFIX, "bwa_prefix", "", "Prefix of BWA index. Should be built from FASTA file provided in 'ref_fasta'");
    args.add_int(Opt::SAMPLE_DIST, "sample_dist", 1000, "Average distance between alignment start locations");
    args.add_int(Opt::RANDOM_SEED, "random_seed", -1, "Seed to use for starting locations. If negative will use system time.");
    args.add_flag(Opt::DYNAMIC_ALIGN, "dynamic_align", "If set will use dynamic programming algorithm to build off previous alignments. Faster for very low sample distances (< 10) on repetitive genomes");
    args.parse_args(argc, argv);

    std::vector< std::vector<u8> > seqs;
    std::ifstream fasta_in(args.get_string(Opt::REF_FASTA));
    std::string fasta_line;
    while (getline(fasta_in, fasta_line)) {
        if (fasta_line[0] == '>') {
            seqs.push_back(std::vector<u8>());
        } else {
            for (char c : fasta_line) {
                seqs.back().push_back(BASE_COMP_B[BASE_BYTES[(u8)c]]);
            }
        }
    }

    BwaFMI fmi(args.get_string(Opt::BWA_PREFIX));

    if (args.get_int(Opt::RANDOM_SEED) < 0) {
        srand(time(NULL));
    } else {
        srand((u32) args.get_int(Opt::RANDOM_SEED));
    }
    
    Timer t;
    for (auto seq : seqs) {
        if (args.get_flag(Opt::DYNAMIC_ALIGN)) {
            dyn_align(fmi, seq, args.get_int(Opt::SAMPLE_DIST));
        } else {
            brute_align(fmi, seq, args.get_int(Opt::SAMPLE_DIST));
        }
    }
    std::cerr << "Self-alignment time: " << t.lap() << "ms" << std::endl;
}
