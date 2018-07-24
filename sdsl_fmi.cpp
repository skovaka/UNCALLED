#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <climits>
#include "sdsl_fmi.hpp"

void str_to_vec(const std::string &str, sdsl::int_vector<8> &vec) {
    vec = sdsl::int_vector<8>(str.size());
    for (size_t i = 0; i < str.size(); i++) {
        vec[i] = BASE_BYTES[(u8) str[i]]+1;
    }
}

SdslFMI::SdslFMI() {
    loaded_ = false;
}

SdslFMI::SdslFMI(const std::string &filename) {
    loaded_ = load_from_file(index_, filename);
}

void SdslFMI::construct(const std::string &seq) {
    sdsl::int_vector<8> vec;
    str_to_vec(seq, vec);
    sdsl::construct_im(index_, vec);
    loaded_ = true;
}

void SdslFMI::save(const std::string &filename) {
    sdsl::store_to_file(index_, filename); // save it
}

Range SdslFMI::get_neighbor(Range r1, u8 base) const {
    Range r2;
    if (r1.start_ <= r1.end_) {
        sdsl::backward_search(index_, r1.start_, r1.end_, base+1, r2.start_, r2.end_);
    }
    return r2;
}

Range SdslFMI::get_full_range(u8 base) const {
    Range r;
    sdsl::backward_search(index_, 0, index_.size()-1, base+1, r.start_, r.end_);
    return r;
}

size_t SdslFMI::sa(size_t i) const {
    return index_[i];
}

size_t SdslFMI::size() const {
    return index_.size();
}

    /*
int main(int argc, char** argv) {
    std::string index_suffix = ".fm9";
    std::string index_file = std::string(argv[1])+index_suffix;
    //sdsl::csa_wt fm_index;
    
    std::cout << argv[1] << std::endl;

    if (!sdsl::load_from_file(fm_index, index_file)) {
        std::ifstream in(argv[1]);

        //For parsing the file
        std::string line, bases;
        
        getline(in, line); //read past header

        while (getline(in, line)) {
            bases += line;
        }

        sdsl::int_vector<8> enc_bases = str_to_vec(bases);
        //u8 *enc_bases = new u8[bases.size()];
        //std::vector<u8> enc_bases(bases.size());
        //for (size_t i = 0; i < bases.size(); i++) {
        //    enc_bases[i] = char_to_base(bases[i])+1;
            //std::cout << (int) enc_bases[i] << std::endl;
        //}

        //std::cout << std::endl << bases << std::endl;
        //fm_index = sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127> >, 512, enc_bases.size()>();
        sdsl::construct_im(fm_index, enc_bases); // generate index
        sdsl::store_to_file(fm_index, index_file); // save it
    }

    sdsl::int_vector<8> pat_vec = str_to_vec(std::string(argv[2]));
    //auto locs = sdsl::locate(fm_index, pat_vec.begin(), pat_vec.end());
    //std::cout << fm_index.size() << " " << locs.size() << std::endl;
    //for (size_t i = 0; i < locs.size(); i++) {
    //    std::cout << locs[i] << std::endl;
    //}

    size_t lo = 0, hi = fm_index.size()-1;
    std::cout << lo << "\t" << hi << std::endl;
    for (size_t i = pat_vec.size()-1; i < pat_vec.size(); i--) {
        sdsl::backward_search(fm_index, lo, hi, pat_vec[i], lo, hi);
        std::cout << lo << "\t" << hi << "\t" << (hi - lo + 1) << std::endl;
        if (lo > hi) {
            break;
        }
    }

    for (size_t i = lo; i <= hi; i++) {
        std::cout << fm_index[i] << std::endl;
    }
    
    cout << "Index c3nstruction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB." << endl;
    cout << "Input search terms and press Ctrl-D to exit." << endl;
    string prompt = "\e[0;32m>\e[0m ";
    cout << prompt;
    string query;
    while (getline(cin, query)) {
        size_t m  = query.size();
        size_t occs = sdsl::count(fm_index, query.begin(), query.end());
        cout << "# of occurrences: " << occs << endl;
        if (occs > 0) {
            cout << "Location and context of first occurrences: " << endl;
            auto locations = locate(fm_index, query.begin(), query.begin()+m);
            sort(locations.begin(), locations.end());
            for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i) {
                cout << setw(8) << locations[i] << ": ";
                if (pre_extract > locations[i]) {
                    pre_extract = locations[i];
                }
                if (locations[i]+m+ post_extract > fm_index.size()) {
                    post_extract = fm_index.size()-locations[i]-m;
                }
                auto s   = extract(fm_index, locations[i]-pre_extract, locations[i]+m+ post_extract-1);
                string pre = s.substr(0, pre_extract);
                s = s.substr(pre_extract);
                if (pre.find_last_of('\n') != string::npos) {
                    pre = pre.substr(pre.find_last_of('\n')+1);
                }
                cout << pre;
                cout << "\e[1;31m";
                cout << s.substr(0, m);
                cout << "\e[0m";
                string context = s.substr(m);
                cout << context.substr(0, context.find_first_of('\n')) << endl;
            }
        }
        cout << prompt;
    }
    cout << endl;

}

    */
