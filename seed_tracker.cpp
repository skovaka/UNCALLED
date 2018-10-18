#include "seed_tracker.hpp"
//#include "base_fmi.hpp"
#include <iostream>
#include <set>

//#define DEBUG

Seed::Seed(u32 read_end, 
           u8 seed_len, 
           float prob, 
           u64 ref_start, 
           u64 ref_end) 
    : read_range_( Range(read_end - seed_len, read_end) ),
      ref_range_(Range(ref_start, ref_end)),
      seed_prob_(prob) {}

void Seed::set_ref_range(u64 end, u8 length) {
    ref_range_.start_ = end - length + 1;
    ref_range_.end_ = end;
}


void Seed::print(std::ostream &out) {
    out << read_range_.start_ << "-" << read_range_.end_ << "\t"
              << ref_range_.start_ << "-" << ref_range_.end_ << "\t"
              << seed_prob_ << "\n";
}

static const ReadAln NULL_ALN = ReadAln();
u8 ReadAln::WIN_LEN;

ReadAln::ReadAln() 
    : evt_st_(1),
      evt_en_(0),
      total_len_(0) {
}

ReadAln::ReadAln(Range ref_st, u32 evt_st)
    : ref_st_(ref_st.start_),
      ref_en_(ref_st),
      evt_st_(evt_st),
      evt_en_(evt_st),
      total_len_(ref_st.length()) {}

ReadAln::ReadAln(const ReadAln &r)
    : ref_st_(r.ref_st_),
      ref_en_(r.ref_en_),
      evt_st_(r.evt_st_),
      evt_en_(r.evt_en_),
      total_len_(r.total_len_) {}


u8 ReadAln::update(ReadAln &new_aln) {
    u8 growth = 0;
    if (new_aln.ref_en_.start_ < ref_en_.end_) {
        if (new_aln.ref_en_.end_ > ref_en_.end_) {
            growth = new_aln.ref_en_.end_ - ref_en_.end_;
            ref_en_ = new_aln.ref_en_;
        } else {
            ref_en_.start_ = new_aln.ref_en_.start_;
        }
    } else {
        growth = new_aln.total_len_;
        ref_en_ = new_aln.ref_en_;
    }

    evt_en_ = new_aln.evt_en_;
    total_len_ += growth;
    return growth;
}

Range ReadAln::ref_range() const {
    return Range(ref_st_, ref_en_.end_);
}

void ReadAln::print(std::ostream &out, bool newline = false, bool print_all = false) const {
    out << total_len_ << "\t";

    out << ref_st_;
    //if (print_all ) {
    //    std::cout << ref_st_.end_ << ":";
    //}

    out << "-" << ref_en_.end_ << "\t" 
               << evt_st_ << "-" 
               << evt_en_;// << " "
               //<< (int) segments_;

    if (newline)
        out << "\n";
}

bool ReadAln::is_valid() {
    return evt_st_ <= evt_en_;
}


bool operator< (const ReadAln &r1, const ReadAln &r2) {
    if (r1.ref_en_.start_ != r2.ref_en_.start_)
        return r1.ref_en_.start_ > r2.ref_en_.start_;

    return r1.evt_en_ > r2.evt_en_;
}

std::ostream &operator<< (std::ostream &out, const ReadAln &a) {
    out << a.ref_st_ << "-" << a.ref_en_.end_ << "\t"
        << a.evt_st_ << "-" << (a.evt_en_ + ReadAln::WIN_LEN) << "\t"
        << a.total_len_;
    return out;
}

SeedTracker::SeedTracker(u64 ref_len, float mean_thresh, float top_thresh, u8 min_aln_len, u8 win_len)
    : ref_len_(ref_len),
      mean_thresh_(mean_thresh),
      top_thresh_(top_thresh),
      min_aln_len_(min_aln_len) {
    ReadAln::WIN_LEN = win_len;
    reset();
}

void SeedTracker::reset() {
    alignments_.clear();
    all_lens_.clear();
    max_len_ = 0;
    len_sum_ = 0;
}

ReadAln SeedTracker::add_seeds(const std::vector<Seed> &seeds) {
    ReadAln top;

    for (size_t i = 0; i < seeds.size(); i++) {
        ReadAln a = add_seed(seeds[i]);

        if (a.total_len_ >= min_aln_len_ && a.total_len_ > max_len_) {


            max_len_ = a.total_len_;
            float mean_len = len_sum_ / alignments_.size();
            float next_len = *std::next(all_lens_.rbegin());

            if ((mean_thresh_ > 0 && max_len_ / mean_len >= mean_thresh_) ||
                (top_thresh_ > 0  && max_len_ / next_len >= top_thresh_)) {
                return a;
            }
        }
    }

    return NULL_ALN;
}

ReadAln SeedTracker::add_seed(Seed r) {
    ReadAln new_aln(r.ref_range_, r.read_range_.start_);

    //Locations sorted by decreasing ref_en_.start
    //Find the largest aln s.t. aln->ref_en_.start <= new_aln.ref_en_.start
    //AKA r1 <= r2
    auto aln = alignments_.lower_bound(new_aln),
         aln_match = alignments_.end();

    u64 e2 = new_aln.evt_en_, //new event aln
        r2 = new_aln.ref_en_.start_; //new ref aln

    while (aln != alignments_.end()) {
        u64 e1 = aln->evt_en_, //old event aln
            r1 = aln->ref_en_.start_; //old ref aln

        //We know r1 <= r2 because of alnation sort order

        bool higher_sup = aln_match == alignments_.end() 
                       || aln_match->total_len_ < aln->total_len_,
             
             in_range = e1 <= e2 && //event aln must increase
                        r2 - r1 <= e2 - e1 && //evt increases more than ref (+ skip)
                        (r2 - r1) >= (e2 - e1) / 12; //evt doesn't increase too much
             
        if (higher_sup && in_range) {
            aln_match = aln;
        } else if (r2 - r1 >= e2) {
            break;
        }

        aln++;
    }

    if (aln_match != alignments_.end()) {
        ReadAln a = *aln_match;

        u32 prev_len = a.total_len_;
        a.update(new_aln);

        if (a.total_len_ != prev_len) {
            len_sum_ += a.total_len_ - prev_len;
            auto l = all_lens_.find(prev_len);
            all_lens_.insert(l, a.total_len_);
            all_lens_.erase(l);
        }

        auto hint = std::next(aln_match);
        alignments_.erase(aln_match);
        alignments_.insert(hint, a);

        #ifdef DEBUG
        new_aln.print(std::cout, true, false);
        #endif

        return a;
    }

    alignments_.insert(new_aln);
    all_lens_.insert(new_aln.total_len_);
    len_sum_ += new_aln.total_len_;

    #ifdef DEBUG
    new_aln.print(std::cout, true, false);
    #endif

    return new_aln;
}

void SeedTracker::print(std::ostream &out, u16 max_out = 10) {
    if (alignments_.empty()) {
        return;
    }

    std::vector<ReadAln> alns_sort(alignments_.begin(),
                                   alignments_.end());

    std::sort(alns_sort.begin(), alns_sort.end(),
              [](const ReadAln &a, const ReadAln &b) -> bool {
                  return a.total_len_ > b.total_len_;
              });

    Range top_ref = alns_sort[0].ref_range();
    double top_len = alns_sort[0].total_len_;

    for (unsigned int i = 0; i < min(max_out, alns_sort.size()); i++) {
        double overlap = top_ref.get_recp_overlap(alns_sort[i].ref_range()),
               len_ratio = top_len / alns_sort[i].total_len_;

        alns_sort[i].print(out, false);
        out << "\t" << len_ratio << "\t" << overlap << "\n";
    }
}

