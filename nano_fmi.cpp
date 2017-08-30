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

#define DEBUG(s)
//#define DEBUG(s) do { std::cerr << s; } while (0)

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

    //Allocate space for other data structures
    bwt_.resize(mer_seq.size());
    mer_f_starts_.resize(model.kmer_count());
    mer_counts_.resize(model.kmer_count());
    mer_tally_.resize(model.kmer_count());

    mer_count_tmp_.resize(model.kmer_count(), 0);

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


float NanoFMI::get_stay_prob(Event e1, Event e2) const {
    double var1 = e1.stdv*e1.stdv, var2 = e2.stdv*e2.stdv;

    double t = (e1.mean - e2.mean) / sqrt(var1/e1.length + var2/e2.length);

    int df = pow(var1/e1.length + var2/e2.length, 2) 
               / (pow(var1/e1.length, 2) / (e1.length-1) 
                    + pow(var2/e2.length, 2) / (e2.length-1));

    boost::math::students_t dist(df);
    double q = boost::math::cdf(boost::math::complement(dist, fabs(t)));

    return q;
}

//Returns the number of occurences of the given k-mer in the BWT up to and
//including the given index
std::list<int> NanoFMI::get_tallies(std::list<mer_id> kmers, int loc) const {
    std::list<int> tallies;
    if (loc < 0)
        return tallies;

    //Closest checkpoint < i
    int cp = (loc / tally_dist_) * tally_dist_; 

    //Check if checkpoint after i is closer
    if (loc - cp > (cp + tally_dist_) - loc 
            && cp + (unsigned) tally_dist_ < bwt_.size())
        cp += tally_dist_;

    int cp_dist = cp - loc; //TODO: just use cp

    for (auto k = kmers.begin(); k != kmers.end(); k++)
        mer_count_tmp_[*k] = mer_tally_[*k][(loc + cp_dist) / tally_dist_];

    if (cp_dist > 0) 
        for (int i = loc+1; i <= loc + cp_dist; i++)
            mer_count_tmp_[bwt_[i]]--;

    else if (cp_dist < 0)
        for (int i = loc; i > loc + cp_dist; i--)
            mer_count_tmp_[bwt_[i]]++;
    
    for (auto k = kmers.begin(); k != kmers.end(); k++) {
        tallies.push_back(mer_count_tmp_[*k]);
    }

    return tallies;
}

std::list<NanoFMI::Range> NanoFMI::get_neigbhors(Range range, std::list<mer_id> kmers) const {
    std::list<Range> results;

    std::list<int> mins = get_tallies(kmers, range.start_ - 1);
    std::list<int> maxs = get_tallies(kmers, range.end_);

    auto kmer = kmers.begin();
    auto min = mins.begin();
    auto max = maxs.begin();
        
    while (kmer != kmers.end()) {
        if (*min < *max) {
            int kmer_st = mer_f_starts_[*kmer];
            results.push_back(Range(kmer_st + *min, kmer_st + *max - 1));
        }

        kmer++;
        min++;
        max++;
    }

    return results;
}

//TODO: Maybe store f as ranges?
NanoFMI::Range NanoFMI::get_full_range(mer_id kmer) const {
    return Range(mer_f_starts_[kmer], mer_f_starts_[kmer] + mer_counts_[kmer] -1 );
}

bool NanoFMI::Range::intersects(const Range &q) const {
    return (start_ >= q.start_ && start_ <= q.end_) ||
           (end_ >= q.start_ && end_ <= q.end_);
}

NanoFMI::Range NanoFMI::Range::split_range(const NanoFMI::Range &r) { 

    Range left;
    if (start_ < r.start_) {
        left = Range(*this);
        left.end_ = r.start_ - 1;
    }

    if (end_ > r.end_) {
        start_ = r.end_ + 1;
    } else {
        start_ = 1;
        end_ = 0;
    }

    return left;
}


NanoFMI::Range& NanoFMI::Range::operator=(const NanoFMI::Range& prev) {
    start_ = prev.start_;
    end_ = prev.end_;
    return *this;
}

NanoFMI::Range::Range(const NanoFMI::Range &prev)
    : start_(prev.start_), 
      end_(prev.end_) {}

NanoFMI::Range::Range() : start_(1), end_(0) {}
NanoFMI::Range::Range(int start, int end) : start_(start), end_(end) {}

bool NanoFMI::Range::same_range(const NanoFMI::Range &q) const {
    return start_ == q.start_ && end_ == q.end_;
}

bool NanoFMI::Range::is_valid() const {
    return start_ <= end_;
}


bool operator< (const NanoFMI::Range &q1, const NanoFMI::Range &q2) {
    if (q1.start_ < q2.start_)
        return true;

    if (q1.start_ == q2.start_ && q1.end_ < q2.end_)
        return true;
    
    return false;
}


//Source constructor
NanoFMI::SeedNode::SeedNode(mer_id kmer, double prob)
    : max_length_(1),
      stay_count_(0),
      kmer_(kmer),
      prob_(prob) {}

//Child constructor
NanoFMI::SeedNode::SeedNode(SeedNode *parent, mer_id kmer, 
                            double prob, bool stay) 
    : max_length_(parent->max_length_ + 1),
      stay_count_(parent->stay_count_ + stay),
      kmer_(kmer),
      prob_(parent->prob_ + prob) {

    parents_.push_front(parent_ptr(parent, stay));
}

//Creates invalid node
NanoFMI::SeedNode::SeedNode()
    : max_length_(0),
      stay_count_(0),
      kmer_(0),
      prob_(0) {}

//Copy constructor
NanoFMI::SeedNode::SeedNode(const SeedNode &s)
    : max_length_(s.max_length_),
      stay_count_(s.stay_count_),
      kmer_(s.kmer_),
      prob_(s.prob_),
      parents_(s.parents_),
      children_(s.children_) {}
    
NanoFMI::SeedNode::~SeedNode() {
    //invalidate();
}
bool NanoFMI::SeedNode::is_valid() {
    return max_length_ > 0;
}

void NanoFMI::SeedNode::invalidate(bool print) {
    if (print) {
        DEBUG("deleting " << max_length_ << " " << parents_.size() << "\n");
    }

    //DEBUG("i\n");
    for (auto p = parents_.begin(); p != parents_.end(); p++) {
        //DEBUG("n\n");
        if (p->first->remove_child(this, true, print)) {
            //DEBUG("v\n");
            delete p->first;
        }
        //DEBUG("a\n");
    }
    //DEBUG("l\n");

    parents_.clear();
    children_.clear();

    max_length_ = stay_count_ = prob_ = 0;
    
    //TODO: remove this node's children?
}

//Returns true if this node should be pruned/erased/deleted/freed
//If should_invalidate = false, will always return false
bool NanoFMI::SeedNode::remove_child(SeedNode *child, bool should_invalidate, bool print) {

    //DEBUG("rm child " << children_.size() << " " << max_length_ << " " << parents_.size() << "\n");

    if (should_invalidate && children_.size() == 1) {
        bool is_source = parents_.empty();
        invalidate(print);

        
        return !is_source;
    }

    //for (auto p = parents_.begin(); p != parents_.end(); p++) {
    //    if (p->first->max_length_ == child->max_length_ + 1) {
    //        if (p->first->remove_child(this, should_invalidate)) {
    //            delete p->first;
    //        }
    //        parents_.erase(p);
    //        break;
    //    }
    //}

    //DEBUG("o\n");
    for (auto c = children_.begin(); c != children_.end(); c++) {
        //DEBUG("v" << children_.size() << "\n");
        //DEBUG(*c << "\n");
        if (*c == child) {
            //DEBUG("e\n");
            children_.erase(c);
            break;
        }
    }
    //DEBUG("d\n");

    //DEBUG("\t" << children_.size()  << "\tchild removed\n");

    return false;
}

NanoFMI::SeedGraph::SeedGraph(const NanoFMI &fmi, 
                              const NormParams &norm_params, 
                              int seed_len, int read_len,
                              double event_prob,
                              double seed_prob,
                              double stay_prob)
    : fmi_(fmi),
      norm_params_(norm_params),
      seed_length_(seed_len),
      cur_event_(read_len), 
      event_prob_(event_prob),
      seed_prob_(seed_prob),
      stay_prob_(stay_prob) {}


std::vector<NanoFMI::Result> NanoFMI::SeedGraph::add_event(Event e) {

    //std::cout << cur_event_ << "\n";
    
    //bool stay = get_stay_prob(events[i], events[i+1]) >= STAY_THRESH;
    //bool stay = false;
    bool stay = true;
    double prob;

    int alph_size = fmi_.model_->kmer_count();
    
    event_kmer_probs_.push_front(new double[alph_size]);
    double *kmer_probs = event_kmer_probs_.front();
    for (int kmer = 0; kmer < alph_size; kmer++)
        kmer_probs[kmer] = fmi_.model_->event_match_prob(e, kmer, norm_params_);
    
    sources_.push_front(std::list<SeedNode *>());

    cur_event_--;

    //std::cout << "Creating children\n";
    //

    //DEBUG("ADDING EVENT " << cur_event_ << " (" << sources_.size() << ")\n");

    int source_count = 0,
        child_count = 0;

    //auto prev_ranges = traversed_ranges.begin(); //Empty on first iteration
    for (auto p = prev_nodes.begin(); p != prev_nodes.end(); p++) {
        //SeedNode *next_node = NULL;

        Range prev_range = p->first;
        SeedNode *prev_node = p->second;


        int prev_kmer = prev_node->kmer_;
        prob = kmer_probs[prev_kmer];

        if (stay && prob >= event_prob_) {
            SeedNode next_node(prev_node, prev_kmer, prob, true);
            SeedNode *nn = add_child(prev_range, next_node);
            child_count++;
        
            //if (nn->parents_.size() > 1) {
            //    //DEBUG("NEW NODE " << nn->max_length_ << " " << nn << " ");
            //    for (auto par = nn->parents_.begin(); par != nn->parents_.end(); par++)
            //        //DEBUG("\t" << par->first->max_length_);
            //    //DEBUG("\n");
            //}
        }

        //DEBUG(prev_node << "(" << prev_node->children_.size() << ")\n");
        
        auto neighbor_itr = fmi_.model_->get_neighbors(prev_kmer);
        std::list<mer_id> next_kmers;

        for (auto n = neighbor_itr.first; n != neighbor_itr.second; n++) 
            if(kmer_probs[*n] >= event_prob_) 
                next_kmers.push_back(*n);


        std::list<Range> next_ranges = fmi_.get_neigbhors(prev_range, next_kmers);
            
        auto next_kmer = next_kmers.begin(); 
        auto next_range = next_ranges.begin();

        
        //std::cout << "Adding nodes\n";
        while (next_kmer != next_kmers.end()) {
            prob = kmer_probs[*next_kmer];

            SeedNode next_node(prev_node, *next_kmer, prob, false);
            SeedNode *nn = add_child(*next_range, next_node);
            child_count++;

            next_kmer++; //Maybe do this in-place up above
            next_range++;
        }
    }

    for (auto p = prev_nodes.begin(); p != prev_nodes.end(); p++) {
        SeedNode *prev_node = p->second;

        if (prev_node->children_.empty()) {
            bool is_source = prev_node->parents_.empty();

            //DEBUG("x" << cur_event_ << "x " << prev_node->max_length_ << "\n");

            if (prev_node->max_length_ == 31) {
                DEBUG("oh no " << cur_event_ << "\n");
            }

            prev_node->invalidate(prev_node->max_length_ == 31);
            
            if (!is_source) {
                delete prev_node;
            } 
        }
    }

    for (mer_id kmer = 0; kmer < alph_size; kmer++) {
        prob = kmer_probs[kmer];
        if (prob >= event_prob_) {
            Range next_range = fmi_.get_full_range(kmer);

            if (next_range.is_valid()) {
                SeedNode next_node(kmer, prob);
                source_count += add_sources(next_range, next_node);
            }
        }
    }
 
    //DEBUG("Added " << source_count << " sources, " << child_count << " children\n");   

    prev_nodes.swap(next_nodes);
    next_nodes.clear();
    
    //print_graph();

    auto r = pop_seeds();


    return r;
}

int NanoFMI::SeedGraph::add_sources(const Range &range, const SeedNode &node) {

    auto lb = next_nodes.lower_bound(range);
    
    auto start = lb;
    while (start != next_nodes.begin() 
           && std::prev(start)->first.intersects(range)){
        start--;
    }

    auto end = lb;
    while (end != next_nodes.end() && end->first.intersects(range)) {
        end++;
    }

    if (start == next_nodes.end() || (start == lb && end == lb)) {

        next_nodes[range] = new SeedNode(node);
        sources_.front().push_back(next_nodes[range]); //TODO: store node, no double query
        return 1;
    }

    std::list<Range> split_ranges;

    Range rr(range); //Right range

    //SeedNode *node_copy;

    for (auto pr = start; pr != end; pr++) {
        Range lr = rr.split_range(pr->first); //Left range

        if (lr.is_valid()) {
            if (lr.same_range(pr->first))
                //DEBUG("WHOOPS B " << pr->second << "\n");
            split_ranges.push_back(lr);
        }

        if (!rr.is_valid())
            break;
    }

    if (rr.is_valid()) 
        split_ranges.push_back(rr);

    for (auto r = split_ranges.begin(); r != split_ranges.end(); r++) {
        next_nodes[*r] = new SeedNode(node);
        sources_.front().push_back(next_nodes[*r]);
    }

    return split_ranges.size();
}

NanoFMI::SeedNode *NanoFMI::SeedGraph::add_child(Range &range, SeedNode &node) {
    
    //Do I need to do this?
    if (!node.is_valid()) {
        return NULL;
    }
    
    //Parent of the node being added
    parent_ptr new_parent = node.parents_.front();

    //Find closest node >= the node being added
    auto lb = next_nodes.lower_bound(range);

    //Node range hasn't been added yet
    if (lb == next_nodes.end() || !lb->first.same_range(range)) {

        //DEBUG("CASE 0\n");
            
        //Allocate memory for a new node
        SeedNode *new_node = new SeedNode(node);


        //DEBUG("ADDING " << new_node << " " << next_nodes.size() << "\n");

        //Store it with it's range
        auto r = next_nodes.insert(lb, std::pair<Range, SeedNode *>(range, new_node));

        //DEBUG("NEW NODE " << new_node << " " << next_nodes.size() << "\n");

        //Update the parent
        new_parent.first->children_.push_front(new_node);

        //Created one node
        return new_node;
    }
    
    //Node associated with same range
    SeedNode *dup_node = lb->second;
    auto dup_parent = dup_node->parents_.begin(); 

    //DEBUG("found dup node of " << dup_node << "\n");

    //This node is for a longer seed
    if (dup_node->max_length_ < node.max_length_) {
        //DEBUG(dup_node << " has new best parent, old still there\n");

        //DEBUG("CASE 1 " << dup_node->max_length_ << " " << node.max_length_ << "\n");

        //if (!dup_node->parents_.empty())
        //    //DEBUG(dup_node->parents_.front().first->max_length_ << " par?\n");

        //Copy seed info
        dup_node->max_length_ = node.max_length_;
        dup_node->stay_count_ = node.stay_count_;
        dup_node->prob_ = node.prob_;

        //Update parents and children
        dup_node->parents_.push_front(new_parent);
        new_parent.first->children_.push_front(dup_node);

        //No new nodes created
        return dup_node;

    } else if (dup_node->max_length_ == node.max_length_) {

        if (dup_node->prob_ < node.prob_) {
            //DEBUG(dup_node << " has new best parent, old is gone\n");

            //If dup_parent has no children in the end, will be invalidated (in add_event)
            dup_parent->first->remove_child(dup_node, false);

            //Erase old parent 
            auto old_loc = dup_node->parents_.erase(dup_parent);

            //Copy seed info
            dup_node->max_length_ = node.max_length_;
            dup_node->stay_count_ = node.stay_count_;
            dup_node->prob_ = node.prob_;

            //Insert new parent in old parent's place
            dup_node->parents_.insert(old_loc, new_parent);
            new_parent.first->children_.push_front(dup_node);
        }

        return dup_node;
    }
    //DEBUG(dup_node << " the hunt is on\n");
   
    //Find where to place the new parent pointer 
    while (dup_parent != dup_node->parents_.end()) {

        //dup_parent is a longer seed
        if (new_parent.first->max_length_ < dup_parent->first->max_length_) {
            dup_parent++;
            continue; //not breaking yet

        //dup_parent is the longest branch shorter than new_parent
        } else if (new_parent.first->max_length_ > dup_parent->first->max_length_)  {
            //DEBUG("CASE 2\n");
            
            //Add parent in the correct location
            dup_node->parents_.insert(dup_parent, new_parent);
            new_parent.first->children_.push_front(dup_node);

        //dup_parent is same length as new_parent

        //and new node has higher probability
        } else if (node.prob_ > dup_node->prob_) {

            //DEBUG("CASE 3\n");
            
            //If dup_parent has no children in the end, will be invalidated (in add_event)
            dup_parent->first->remove_child(dup_node, false);

            //Erase old parent 
            auto old_loc = dup_node->parents_.erase(dup_parent);

            //Insert new parent in old parent's place
            dup_node->parents_.insert(old_loc, new_parent);
            new_parent.first->children_.push_front(dup_node);

        } //If new node has lower prob, new node isn't stored in any way

        //Breaking before dup_parent reaches end, so next if statement won't execute
        break;
    }

    //Smallest parent, put at the end
    if (dup_parent == dup_node->parents_.end()) {
        //DEBUG("CASE 4\n");    
        dup_node->parents_.push_back(new_parent);
        new_parent.first->children_.push_front(dup_node);
    }

    return dup_node;
}



std::vector<NanoFMI::Result> NanoFMI::SeedGraph::pop_seeds() { //Big result gathering loop

    //Is there a long enough seed?
    if (sources_.size() < seed_length_)
        return std::vector<NanoFMI::Result>();

    //DEBUG("SEED POP\n");

    //TODO: Store prev_nodes like this?
    std::unordered_map<SeedNode *, Range> node_ranges;
    for (auto n = prev_nodes.begin(); n != prev_nodes.end(); n++)
        node_ranges[n->second] = n->first;


    std::list<SeedNode *> &aln_ends = sources_.back(),
                          &next_sources = *std::next(sources_.rbegin());

    while (!aln_ends.empty()) {
        SeedNode *aln_en = aln_ends.front();

        if (!aln_en->is_valid()) {
            aln_ends.pop_front();
            delete aln_en;
            continue;
        }

        std::list<SeedNode *> to_visit;
        for (auto c = aln_en->children_.begin(); c != aln_en->children_.end(); c++) {
            (*c)->parents_.clear();
            next_sources.push_back(*c);
            to_visit.push_back(*c);
        }
        
        double prev_len = aln_en->max_length_;
        auto kmer_probs = event_kmer_probs_.rbegin();

        while (!to_visit.empty()) {


            SeedNode *n = to_visit.front();
            to_visit.pop_front();

            //DEBUG("updating " << n->max_length_);

            if (n->max_length_ != prev_len) {
                prev_len = n->max_length_;
                kmer_probs++;
            }

            //Check if parents are compatible
            if (n->parents_.size() > 1) {
                auto seed_parent = n->parents_.begin(), //where we came from
                     next_parent = std::next(seed_parent); //seed with closest length

                //Parents now have the same length
                if (seed_parent->first->max_length_ == next_parent->first->max_length_) {
                    SeedNode *to_erase;

                    //Erase where we came from
                    if (seed_parent->first->prob_ < next_parent->first->prob_) {
                        to_erase = seed_parent->first; 
                        n->parents_.erase(seed_parent);

                    //Erase other seed branch
                    } else {
                        to_erase = next_parent->first;
                        n->parents_.erase(next_parent);
                    }

                    //Delete branch (if it's not a source)
                    if (to_erase->remove_child(n, true)) {
                        delete to_erase;
                    }
                }
            }

            if (n->max_length_ == seed_length_) {
                if (n->stay_count_ < 12) {
                    Range r = node_ranges[n];
                    for (int s = r.start_; s <= r.end_; s++) {
                        Result r;
                        r.qry_start = cur_event_;// seed_end - total_length + 1;
                        r.qry_end = cur_event_ + n->max_length_ - 1;
                        r.ref_start = fmi_.suffix_ar_[s];
                        r.ref_end = fmi_.suffix_ar_[s] + n->max_length_ - n->stay_count_ - 1;
                        r.prob = n->prob_ / n->max_length_;
                        //results.push_back(r);

                        std::cout << "= rev\t" 
                                  << r.qry_start << "-" << r.qry_end << "\t"
                                  << r.ref_start << "-" << r.ref_end << "\t"
                                  << r.prob << " =\n";
                    }
                }
            } else {
                for (auto c = n->children_.begin(); c != n->children_.end(); c++) {
                    to_visit.push_back(*c);
                }
            }

            n->max_length_--;
            
            if (!n->parents_.empty()) {
                parent_ptr new_parent = n->parents_.front();

                n->prob_ = new_parent.first->prob_ + (*kmer_probs)[n->kmer_];
                n->stay_count_ = new_parent.first->stay_count_ + new_parent.second;
            } else {
                n->prob_ = (*kmer_probs)[n->kmer_];
                n->stay_count_ = 0;
            }
        }
        
        aln_ends.pop_front();
        delete aln_en;
    }

    sources_.pop_back();

    return std::vector<NanoFMI::Result>();

}

//work on this
void NanoFMI::SeedGraph::print_graph() {
    std::set< std::pair<int, SeedNode *> > to_visit;
    auto event_sources = sources_.rbegin();

    std::cout << "== nodeptr \"event\" max_len #parents children... (" 
              << sources_.size() << " events) ==\n";

    for (auto s = event_sources->begin(); s != event_sources->end(); s++) 
        to_visit.insert( std::pair<int, SeedNode *> (0, *s) );

    bool iter_sources = true;

    int source_count = 0;

    while (!to_visit.empty()) {

        if (iter_sources) {
            source_count++;
            event_sources++;
        }

        if (event_sources != sources_.rend()) {
            for (auto s = event_sources->begin(); s != event_sources->end(); s++) 
                to_visit.insert(std::pair<int, SeedNode *>(source_count, *s));
        } else {
            iter_sources = false;
        }

        auto p = to_visit.begin();
        int event = p->first;
        SeedNode *n = p->second;
        to_visit.erase(p);
        
        std::cout << n << " " << event << " " 
                  << n->max_length_ << " " 
                  << n->parents_.size();

        for (auto c = n->children_.begin(); c != n->children_.end(); c++) {
            std::cout << " " << *c;
            to_visit.insert(std::pair<int, SeedNode *>(event + 1, *c));
        }

        std::cout << "\n";

    }

    std::cout << "== end ==\n";
}
