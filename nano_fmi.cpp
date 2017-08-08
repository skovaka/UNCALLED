#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stddef.h>
#include <list>
#include <set>
#include <utility>
#include <unordered_map>
#include "timer.h"
#include "nano_fmi.hpp"
#include "boost/math/distributions/students_t.hpp"

//Reads a model directly from a file and creates the FM index from the given reference
NanoFMI::NanoFMI(KmerModel &model, std::vector<mer_id> &mer_seq, int tally_dist) {

    model_ = &model;
    mer_seq_ = &mer_seq;
    tally_dist_ = tally_dist;

    //For outputting time
    Timer timer;

    //Init suffix array
    //Not using suffix_ar instance var speeds up sorting significantly
    std::vector<unsigned int> suffix_ar(mer_seq.size());
    for (unsigned int i = 0; i < suffix_ar.size(); i++) 
        suffix_ar[i] = i;

    std::cerr << "SA init time: " << timer.lap() << "\n";

    //Create the suffix array
    std::sort(suffix_ar.begin(), suffix_ar.end(), *this);
    suffix_ar_.swap(suffix_ar);

    std::cerr << "SA sort time: " << timer.lap() << "\n";

    //Allocate space for other datastructures
    bwt_.resize(mer_seq.size());
    mer_f_starts_.resize(model.kmer_count());
    mer_counts_.resize(model.kmer_count());
    mer_tally_.resize(model.kmer_count());

    for (mer_id i = 0; i < model.kmer_count(); i++)
        mer_tally_[i].resize((mer_seq.size() / tally_dist_) + 1, -1);
    
    std::cerr << "FM init time: " << timer.lap() << "\n";

    int tally_mod = tally_dist_;
    
    //Single pass to generate BWT and other datastructures
    for (unsigned int i = 0; i < suffix_ar_.size(); i++) {
        
        //Fill in BWT
        if (suffix_ar_[i] > 0)
            bwt_[i] = mer_seq[suffix_ar_[i]-1];
        else
            bwt_[i] = mer_seq[suffix_ar_[suffix_ar_.size()-1]];

        //Update 6-mer counts
        mer_counts_[bwt_[i]]++;
        
        //Update tally array
        if (tally_mod == tally_dist_) {
            for (mer_id j = 0; j < model.kmer_count(); j++)
                mer_tally_[j][i / tally_dist_] = mer_counts_[j];
            tally_mod = 0;
        }
        tally_mod += 1;
    }

    std::cerr << "FM build time: " << timer.lap() << "\n";
    
    //Compute start locations for F array
    for (mer_id i = 0; i < model.kmer_count(); i++) {
        mer_f_starts_[i] = 1;
        for (mer_id j = 0; j < model.kmer_count(); j++)
            if (model_->compare_kmers(i, j) > 0)
                mer_f_starts_[i] += mer_counts_[j];
    }
    
    //Fill in last entry in tally array if needed
    if (mer_seq.size() % tally_dist_ == 0)
        for (mer_id i = 0; i < model.kmer_count(); i++)
            mer_tally_[i][mer_tally_[i].size()-1] = mer_counts_[i];

}

//Returns true if the suffix of *mer_seq_tmp starting at rot1 is less than that
//starting at rot2. Used to build suffix array.
bool NanoFMI::operator() (unsigned int rot1, unsigned int rot2) {
    
    int c;
    for (unsigned int i = 0; i < mer_seq_->size(); i++) {
        c = model_->compare_kmers(mer_seq_->at(rot1+i), mer_seq_->at(rot2+i));
        
        if (c == 0)
            continue;

        if (c < 0)
            return true;

       return false;
    }

    return false;
}


float NanoFMI::get_stay_prob(Event e1, Event e2) {
    double var1 = e1.stdv*e1.stdv, var2 = e2.stdv*e2.stdv;

    double t = (e1.mean - e2.mean) / sqrt(var1/e1.length + var2/e2.length);

    int df = pow(var1/e1.length + var2/e2.length, 2) 
               / (pow(var1/e1.length, 2) / (e1.length-1) 
                    + pow(var2/e2.length, 2) / (e2.length-1));

    boost::math::students_t dist(df);
    double q = boost::math::cdf(boost::math::complement(dist, fabs(t)));

    return q;
}

//Returns the distance from a BWT index to the nearest tally array checkpoint
int NanoFMI::tally_cp_dist(int i) {
    int cp = (i / tally_dist_)*tally_dist_; //Closest checkpoint < i

    //Check if checkpoint after i is closer
    if (i - cp > (cp + tally_dist_) - i && cp + (unsigned) tally_dist_ < bwt_.size())
        cp += tally_dist_;

    return cp - i;
}

//Returns the number of occurences of the given k-mer in the BWT up to and
//including the given index
int NanoFMI::get_tally(mer_id c, int i) {
    if (i < 0)
        return -1;
    int cp_dist = tally_cp_dist(i);
    int tally = mer_tally_[c][(i + cp_dist) / tally_dist_];

    if (cp_dist > 0) {
        for (int j = i+1; j <= i + cp_dist; j++)
            if (bwt_[j] == c)
                tally--;

    } else if (cp_dist < 0) {
        for (int j = i; j > i + cp_dist; j--)
            if (bwt_[j] == c)
                tally++;
    }

    return tally;
}

//Aligns the vector of events to the reference using LF mapping
//Returns the number of exact alignments
std::vector<NanoFMI::Result> NanoFMI::lf_map(
                    std::vector<Event> &events, int map_start, 
                    int seed_len, NormParams norm_params) {

    //float EVENT_THRESH = 0.00025,
    //      SEED_THRESH  = 0.05,
    //      STAY_THRESH  = 0.01;

    float EVENT_THRESH = -9.2103,
          SEED_THRESH = -3.75, //-5.298,
          STAY_THRESH = -5.298;


    //Match the first event

    double event_probs[map_start + 1][model_->kmer_count()];

    std::cout << "Getting probs\n";

    for (int e = map_start; e >= 0; e--) {
        for (mer_id i = 0; i < model_->kmer_count(); i++) {
            event_probs[e][i] = model_->event_match_prob(events[e], i, norm_params);
        }
    }

    std::cout << "Got probs\n";

    std::list< std::set<Query> > queries(1);

    std::list< std::list<mer_id> > matched_mers(1);
    std::set<Query> finished;

    //double prob;

    for (int seed_end = map_start; seed_end >= 2226; seed_end--) {


        for (mer_id i = 0; i < model_->kmer_count(); i++) {
            if (event_probs[map_start][i] >= EVENT_THRESH) {
                Query q(*this, i, event_probs[map_start][i]);
                update_queries(queries.back(), q);
            }
        }

        bool mer_matched = true,
             stay = false;

        auto prev_queries = queries.rbegin();
        std::cout << "Starting loop: " << seed_end << " " << queries.size() << "\n";
        

        //Match all other events backwards
        int i = seed_end - 1;
        while (i >= 0 && mer_matched) {

            if (&(*prev_queries) == &(queries.front())) {
                queries.push_front(std::set<Query>());
            }
            
            auto next_queries = prev_queries;
            next_queries++;

            mer_matched = false;

            //stay = get_stay_prob(events[i], events[i+1]) >= STAY_THRESH;
            stay = true;

            for (auto pq = prev_queries->begin(); pq != prev_queries->end(); pq++) {

                if (!pq->children_->empty())
                    continue;

                if (stay && event_probs[i][pq->k_id_] >= EVENT_THRESH) {
                    Query nq(*pq, event_probs[i][pq->k_id_]);
                    update_queries(*next_queries, nq);
                }
                
                auto neighbors = model_->get_neighbors(pq->k_id_);

                for (auto n = neighbors.first; n != neighbors.second; n++) {

                    if(event_probs[i][*n] >= EVENT_THRESH) {

                        Query nq(*pq, *n, event_probs[i][*n]);

                        if (nq.is_valid()) {
                            if ( nq.match_len_ < seed_len ) {

                                update_queries(*next_queries, nq);

                                mer_matched = true;

                            } else {
                                finished.insert(nq);
                            }
                        }
                    }
                }
            }


            prev_queries++;

            //Stop search if not matches found
            if (!mer_matched)
                break;

            i--;
        }

        break;

    }

    int event_count = 0;
    for (auto s = queries.begin(); s != queries.end(); s++) {
        //for (auto q = s->begin(); q != s->end(); q++) {
        //    std::cout << q->size() << "\n";
        //}
        //std::cout << s->size() << "\n";
        event_count += s->size();
    }
    std::cout << event_count << " number of events\n";
    std::cout << "fin: " << finished.size() << "\n";
    


    std::vector<Result> results;
    for (auto qry = finished.begin(); qry != finished.end(); qry++) {
        for (int s = qry->start_; s <= qry->end_; s++) {
            
            const Query *q = &(*qry);
            int length = 2;
            while (!q->parents_->empty()) {
                q = q->parents_->front();
                std::cout << q->parents_->size() << "\n";
                length++;
            }

            int ref_start = suffix_ar_[s];
            int ref_end = suffix_ar_[s] + qry->match_len_ - 1;
            std::cout << ref_start << "-" << ref_end << " " 
                      << qry->parents_->size() << " " << qry->prob_/length << "\n";


        }
    }

    return results;
}

int NanoFMI::update_queries(std::set<NanoFMI::Query> &queries, NanoFMI::Query &new_qry) {

    if (queries.empty()) {
        queries.insert(new_qry);
        return 1;
    }

    auto lb = queries.lower_bound(new_qry);

    if (lb == queries.end()) {
        queries.insert(lb, new_qry);
        return 1;
    }

    if (lb->same_range(new_qry)) {
        if (new_qry.prob_ > lb->prob_) {
            queries.erase(lb);
            queries.insert(new_qry);
        }
        //lb->parents_->insert(lb->parents_->end(), 
        //                     new_qry.parents_->begin(), 
        //                     new_qry.parents_->end());
        return 0;        
    }

    auto end = lb;
    while (end != queries.end() && end->intersects(new_qry)) {
        end++;
    }

    auto start = lb;
    while (start != queries.begin() && std::prev(start)->intersects(new_qry)){
        start--;
    }
    

    if (start == queries.end() || (start == lb && end == lb)) {
        queries.insert(lb, new_qry);
        return 1;
    }

    std::list<Query> split_queries;
    for (auto sq = start; sq != end; sq++) {
        Query sq2(*sq);
        split_queries.push_back(sq2.split_query(new_qry)); //q is right query now
    }
    split_queries.push_back(new_qry);


    int new_count = 0;
    for (auto sq = split_queries.begin(); sq != split_queries.end(); sq++) {
        if (sq->is_valid()) {
            queries.insert(lb, *sq);
            new_count++;
        }
    }
    
    return new_count;
}

bool NanoFMI::Query::intersects(const Query &q) const {
    return (start_ >= q.start_ && start_ <= q.end_) ||
           (end_ >= q.start_ && end_ <= q.end_);
}

NanoFMI::Query NanoFMI::Query::split_query(NanoFMI::Query &q) {
    Query left(fmi_);

    if (start_ > q.start_) {
        left = Query(q);
        left.end_ = start_ - 1;
    }

    if (end_ < q.end_) {
        q.start_ = end_ + 1;
    }

//    parents_->insert(parents_->begin(), q.parents_->begin(), q.parents_->end());

    return left;
}


NanoFMI::Query& NanoFMI::Query::operator=(const NanoFMI::Query& prev) {
    fmi_ = prev.fmi_;
    k_id_ = prev.k_id_;
    start_ = prev.start_;
    end_ = prev.end_;
    prob_ = prev.prob_;
    match_len_ = prev.match_len_;
    parents_ = new std::list<const Query *> (*prev.parents_);
    children_ = new std::list<const Query *> (*prev.children_);
    return *this;
}

NanoFMI::Query::Query(NanoFMI &fmi)
    : fmi_(fmi), 
      k_id_(0),
      start_(0), 
      end_(0), 
      match_len_(0),
      prob_(0),
      parents_(NULL),
      children_(NULL) {}

//Initial match constructor
NanoFMI::Query::Query(NanoFMI &fmi, mer_id k_id, double prob)
    : fmi_(fmi), 
      k_id_(k_id),
      start_(fmi.mer_f_starts_[k_id]), 
      end_(fmi.mer_f_starts_[k_id] + fmi.mer_counts_[k_id] - 1), 
      match_len_(1),
      prob_(prob),
      parents_(new std::list<const Query *>),
      children_(new std::list<const Query *>) {}

//"next" constructor
NanoFMI::Query::Query(const NanoFMI::Query &prev, mer_id k_id, double prob)
    : fmi_(prev.fmi_),
      k_id_(k_id),
      parents_(new std::list<const Query *>),
      children_(new std::list<const Query *>) {

    int min = fmi_.get_tally(k_id, prev.start_ - 1),
        max = fmi_.get_tally(k_id, prev.end_);

    parents_->push_back(&prev);
    prev.children_->push_back(this);
    
    //Candidate k-mer occurs in the range
    if (min < max) {
        start_ = fmi_.mer_f_starts_[k_id] + min;
        end_ = fmi_.mer_f_starts_[k_id] + max - 1;
        match_len_ = prev.match_len_ + 1;
        prob_ = prev.prob_ + prob;
    } else {
        start_ = end_ = 0;

    }
}

//"stay" constructor
NanoFMI::Query::Query(const NanoFMI::Query &prev, double prob)
    : fmi_(prev.fmi_), 
      k_id_(prev.k_id_),
      start_(prev.start_), 
      end_(prev.end_), 
      match_len_(prev.match_len_),
      prob_(prev.prob_ + prob),
      parents_(new std::list<const Query *>),
      children_(new std::list<const Query *>) {
    
    parents_->push_back(&prev);
    prev.children_->push_back(this);
        
}

//Copy constructor
NanoFMI::Query::Query(const NanoFMI::Query &prev)
    : fmi_(prev.fmi_), 
      k_id_(prev.k_id_),
      start_(prev.start_), 
      end_(prev.end_), 
      match_len_(prev.match_len_),
      prob_(prev.prob_) {
    
    if (prev.parents_ != NULL)
        parents_ = new std::list<const Query *>(*prev.parents_);
    else
        parents_ = NULL;

    if (prev.children_ != NULL)
        children_ = new std::list<const Query *>(*prev.children_);
    else
        children_ = NULL;

} 

NanoFMI::Query::~Query() {
    //delete parents_;
    //delete children_;
}

bool NanoFMI::Query::add_results(std::vector<NanoFMI::Result> &results, 
                                 int query_end, double min_prob) const {
    //if (avg_prob() >= min_prob) {
    //    Result r;
    //    for (int s = start_; s <= end_; s++) {
    //        r.qry_start = query_end - match_len_ - stays_ + 1;
    //        r.qry_end = query_end;
    //        r.ref_start = fmi_.suffix_ar_[s];
    //        r.ref_end = fmi_.suffix_ar_[s] + match_len_ - 1;
    //        r.prob = avg_prob();
    //        results.push_back(r);
    //    }
    //    return true;
    //}
    //return false;
    return false;
}

bool NanoFMI::Query::same_range(const NanoFMI::Query &q) const {
    return start_ == q.start_ && end_ == q.end_;
}

bool NanoFMI::Query::is_valid() const {
    return start_ <= end_ && (start_ != 0 || end_ != 0);
}


void NanoFMI::Query::print_info() const {
}


bool operator< (const NanoFMI::Query &q1, const NanoFMI::Query &q2) {
    if (q1.start_ < q2.start_)
        return true;

    if (q1.start_ == q2.start_ && q1.end_ < q2.end_)
        return true;
    
    return false;
}
