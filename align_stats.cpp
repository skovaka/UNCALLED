#include <iostream>
#include <string>
#include <list>
#include "base_fmi.hpp"
#include "sdsl_fmi.hpp"
#include "bwa_fmi.hpp"
#include "fmi.hpp"
#include "timer.hpp"

void brute_align(FMI &fmi, const std::vector<u8> &bases) {
    for (size_t j = bases.size()-1; j < bases.size(); j--) {
        Range r = fmi.get_full_range(bases[j]);
        size_t i = j-1;
        for (; i < bases.size() && r.length() > 1; i--) {
            //std::cout << r.length() << "\t";
            r = fmi.get_neighbor(r, bases[i]);
        }
        //std::cout << r.length() << "\n";
        //std::cout << (j - i) << "\n";
    }
}

void dyn_align(FMI &fmi, const std::vector<u8> &bases) {//, std::ofstream &outfile) {
    std::list<Range> ranges;

    u64 rangelen;

    //Align from first base until fm range reaches 1
    ranges.push_back(fmi.get_full_range(bases[0]));
    for (u64 i = 1; i < bases.size() && ranges.back().length() > 1; i++) {
        rangelen = ranges.back().length();//length of previous range
        ranges.push_back( fmi.get_neighbor(ranges.back(), bases[i]) );

        //Output stuff
        //outfile.write((char *) &rangelen, sizeof(unsigned int));
        std::cout << rangelen << "\t";
    }
    
    rangelen = ranges.back().length();

    //Output
    //outfile.write((char *) &rangelen, sizeof(unsigned int));
    std::cout << rangelen << std::endl;
    //std::cout << ranges.size() << "\n";

    //Find alignments starting from each base
    for (u64 i = 1; i < bases.size(); i++) {
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
            std::cout << m->length() << "\t";

            j++;
        }

        //If m != end, right+left must have reached 0 before end
        //Means removing first base didn't add new alignments

        //Output previously found ranges again
        for (; m != ranges.end(); m++) {
            std::cout << m->length() << "\t";
        }

        //Don't output previous ranges again, just mark it
        //if (m != ranges.end()) {
        //    std::cout << "x\t";
        //}

        j = i + ranges.size();

        //bool added = false;
        //
        //Removing first base DID add new alignments
        while (ranges.back().length() > 1 && j < bases.size()) {
            r = fmi.get_neighbor(ranges.back(), bases[j]);
            if (!r.is_valid()) { 
                break;
            }
            ranges.push_back(r);
            std::cout << r.length() << "\t";
            j++;
            //added = true;

            //Output
        }

        std::cout << std::endl;

        //if (added) {
            //for (auto k = ranges.begin(); k != ranges.end(); k++) {
             //   std::cout << k->length() << "\t";
                //rangelen = (unsigned int) k->length();
                //outfile.write((char *) &rangelen, sizeof(unsigned int));
            //}
            //std::cout << "\n";

        //i_mod += 1;
        //if (i_mod == flush_checkpoint) {
        //    i_mod = 0;
        //    outfile.flush();
        //    std::cout.flush();
        //}
        //}

        //std::cout << "\n";
        //std::cout << ranges.size() << "\n";
    }
}

int main(int argc, char** argv) {
    std::string ref_fname(argv[1]), 
                index_fname(argv[2]),
                //out_fname(argv[3]), 
                fwd_str, rev_str;

    std::cerr << "Reading reference\n";
    std::ifstream ref_file(ref_fname);
    parse_fasta(ref_file, fwd_str, rev_str, false);
    std::vector<u8> bases = seq_to_bases(fwd_str);
    for (u64 i = 0; i < bases.size(); i++) {
        bases[i] = BASE_COMP_B[bases[i]];
    }
    rev_str.erase();
    fwd_str.erase();

    std::cerr << "Reading index\n";
    BwaFMI fmi(index_fname);
    //SdslFMI fmi(index_fname);
    //std::ifstream index_file(index_fname);
    //BaseFMI fmi(index_file, 128);

    //std::ofstream outfile(out_fname, std::ios::out | std::ios::binary);

    std::cerr << "Aligning\n";
    
    Timer t;
    dyn_align(fmi, bases);//, outfile);
    //outfile.close();
    std::cerr << t.lap() << "\n";
}
