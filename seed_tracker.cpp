#include "seed_tracker.hpp"
//#include "base_fmi.hpp"
#include <iostream>
#include <set>

//#define DEBUG

ReadAln::ReadAln() 
    : evt_st_(-1),
      total_len_(0) {
}

ReadAln::ReadAln(Range ref_st, unsigned int evt_st, float prob)
    : ref_en_(ref_st),
      evt_st_(evt_st),
      ref_st_(ref_st.end_),
      evt_en_(evt_st),
      total_len_(ref_st.length()),
      prob_sum_(prob),
      segments_(1) {}

ReadAln::ReadAln(const ReadAln &r)
    : ref_en_(r.ref_en_),
      evt_st_(r.evt_st_),
      ref_st_(r.ref_st_),
      evt_en_(r.evt_en_),
      total_len_(r.total_len_),
      prob_sum_(r.prob_sum_),
      segments_(r.segments_) {}


void ReadAln::update_next(ReadAln &new_loc) const {
    new_loc.ref_st_ = ref_st_;
    new_loc.evt_st_ = evt_st_;
    new_loc.prob_sum_ += prob_sum_;
    //new_loc.segments_ = segments_ + 1;

    if (ref_en_.intersects(new_loc.ref_en_)) {

        if (new_loc.ref_en_.end_ > ref_en_.end_) {
            new_loc.total_len_ 
                = total_len_ + 
                  new_loc.ref_en_.end_ - 
                  ref_en_.end_;
        } else {
            new_loc.total_len_ = total_len_;
            new_loc.ref_en_.end_ = ref_en_.end_;
        }

        new_loc.segments_ = segments_;

    } else {
        new_loc.total_len_ += total_len_;
        new_loc.segments_ = segments_ + 1;
    }

}

Range ReadAln::ref_range() const {
    return Range(ref_st_, ref_en_.end_);
}

void ReadAln::print(std::ostream &out, bool newline = false, bool print_all = false) const {
    out << total_len_ << " ";

    out << ref_st_;
    //if (print_all ) {
    //    std::cout << ref_st_.end_ << ":";
    //}

    out << "-" << ref_en_.end_ << " " 
               << evt_st_ << "-" 
               << evt_en_;// << " "
               //<< (int) segments_;

    if (newline)
        out << "\n";
}

bool ReadAln::is_valid() {
    return evt_st_ >= 0;
}

bool operator< (const ReadAln &r1, const ReadAln &r2) {
    if (r1.ref_en_.start_ != r2.ref_en_.start_)
        return r1.ref_en_.start_ > r2.ref_en_.start_;

    return r1.evt_en_ > r2.evt_en_;
}

SeedTracker::SeedTracker(unsigned int read_length) {
    read_length_ = read_length;
    longest_seed = NULL;
}

void SeedTracker::reset() {
    locations.clear();
    longest_seed = NULL;
}

ReadAln SeedTracker::add_seeds(const std::vector<Result> &seeds) {
    ReadAln top;

    for (size_t i = 0; i < seeds.size(); i++) {
        ReadAln a = add_seed(seeds[i]);

        if (!top.is_valid() || top.total_len_ < a.total_len_) {
            top = a;
        }
    }

    return top;
}

ReadAln SeedTracker::add_seed(Result r) {
    ReadAln new_loc(r.ref_range_, r.read_range_.start_, r.seed_prob_);
    //ReadAln new_loc(r.ref_range_, 
    //                read_length_ - r.read_range_.start_ - 1, 
    //                r.seed_prob_);

    auto loc = locations.lower_bound(new_loc),
         loc_match = locations.end();

    unsigned int aln_length = read_length_ - new_loc.evt_st_;
    //std::cout << aln_length << std::endl;

    while (loc != locations.end()) {
        bool higher_sup = loc_match == locations.end() 
                       || loc_match->total_len_ < loc->total_len_,
             
             in_range = loc->ref_en_.start_ + new_loc.evt_en_
                         >= new_loc.ref_en_.start_ + loc->evt_en_;
        
        //ASSERT loc->ref_st_.end_ >= new_loc.ref_st_.end_

        //if (loc->ref_st_.end_ < new_loc.ref_st_.end_) {
        //    std::cout << "WEIRD" << std::endl;
        //}

        if (higher_sup && in_range) {
            loc_match = loc;
            //std::cout << (loc->ref_st_.end_ - new_loc.ref_st_.end_) << "\t" 
            //          << (loc->evt_st_ - new_loc.evt_st_) << std::endl;
        //} else if (loc->ref_st_.end_ - new_loc.ref_st_.end_ >= aln_length) {
        //    break;
        }

        loc++;
    }

    if (loc_match != locations.end()) {

        loc_match->update_next(new_loc);
        locations.erase(loc_match);
        locations.insert(new_loc);


        #ifdef DEBUG
        new_loc.print(std::cout, true, false);
        #endif

    } else if (loc_match == locations.end()) {
        locations.insert(new_loc);

        #ifdef DEBUG
        new_loc.print(std::cout, true, false);
        #endif
    }

    if (longest_seed == NULL || new_loc.total_len_ > longest_seed->total_len_) {
        longest_seed = &new_loc;
    }

    //std::cerr << new_loc.ref_en_ << " " << r.ref_range_.length() << "\n";

    return new_loc;
}

std::vector<ReadAln> SeedTracker::get_alignments(unsigned int min_len = 1) {
    std::map<unsigned int, ReadAln> sorted_alns;

    for (auto a = locations.begin(); a != locations.end(); a++) {
        if (a->total_len_ >= min_len) {
            sorted_alns[a->total_len_] = *a;
        }
    }

    std::vector<ReadAln> ret;
    for (auto a = sorted_alns.rbegin(); a != sorted_alns.rend(); a++) {
        ret.push_back(a->second);
    }

    return ret;
}

bool SeedTracker::check_ratio(const ReadAln &aln, double ratio) {

    Range top_ref = aln.ref_range();
    unsigned int max_len = (aln.total_len_) / ratio;

    for (auto a = locations.begin(); a != locations.end(); a++) {
        if (top_ref.get_recp_overlap(a->ref_range()) == 0 && 
            a->total_len_ >= max_len) {
            
            return false;
        }
    }

    return true;
}

void SeedTracker::print(std::ostream &out, std::string strand, size_t max_out = 10) {

    std::vector<ReadAln> alns = get_alignments(1);

    //std::cout << strand << " " << top_ratio(alns) << "\n";

    if (alns.empty()) {
        return;
    }

    Range top_ref = alns[0].ref_range();
    double top_len = alns[0].total_len_;

    for (unsigned int i = 0; i < min(max_out, alns.size()); i++) {
        double overlap = top_ref.get_recp_overlap(alns[i].ref_range()),
               len_ratio = top_len / alns[i].total_len_;

        out << strand << " ";
        alns[i].print(out, false);
        out << "\t" << len_ratio << "\t" << overlap << "\n";
    }
}

