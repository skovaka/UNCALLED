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
#include "seed_graph.hpp"
#include "kmer_fmi.hpp"
#include "boost/math/distributions/students_t.hpp"

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)


//Source constructor
SeedGraph::Node::Node(mer_id kmer, double prob)
    : length_(1),
      stay_count_(0),
      skip_count_(0),
      ignore_count_(0),
      kmer_(kmer),
      event_prob_(prob),
      seed_prob_(prob)

      #ifdef DEBUG_PROB
      , min_evt_prob_(prob)
      #endif
      
       {}

//Child constructor
SeedGraph::Node::Node(Node *parent, mer_id kmer, 
                      double prob, Type type) 
    : length_(parent->length_ + 1),
      stay_count_(parent->stay_count_),
      skip_count_(parent->skip_count_),
      ignore_count_(parent->ignore_count_),
      kmer_(kmer),
      event_prob_(prob),
      seed_prob_(parent->seed_prob_ + prob) {

    parents_.push_front(parent_ptr(parent, type));
    update_info();
}

//Creates invalid node
SeedGraph::Node::Node()
    : length_(0),
      stay_count_(0),
      skip_count_(0),
      ignore_count_(0),
      kmer_(0),
      event_prob_(0),
      seed_prob_(0) {}

//Copy constructor
SeedGraph::Node::Node(const Node &s)
    : length_(s.length_),
      stay_count_(s.stay_count_),
      skip_count_(s.stay_count_),
      ignore_count_(s.ignore_count_),
      kmer_(s.kmer_),
      event_prob_(s.event_prob_),
      seed_prob_(s.seed_prob_),

      #ifdef DEBUG_PROB
      min_evt_prob_(s.min_evt_prob_),
      #endif

      parents_(s.parents_),
      children_(s.children_) {}
    
bool SeedGraph::Node::is_valid() {
    return length_ > 0;
}

int SeedGraph::Node::seed_len() {
    return length_;
}

int SeedGraph::Node::match_len() {
    return seed_len() - stay_count_ - ignore_count_ + skip_count_;
}

//TODO: Really not a very good name
bool SeedGraph::Node::better_than(const Node *node) {
    return node->ignore_count_ > ignore_count_ ||
            (node->ignore_count_ == ignore_count_ && 
             node->seed_prob_ < seed_prob_);
}

bool SeedGraph::Node::should_report(const AlnParams &params) {
    return parents_.front().second == Node::Type::MATCH && 
           stay_count_ <= params.max_stay_frac_ * length_ && 
           seed_prob_ >= params.min_seed_pr_ * length_;
}

void SeedGraph::Node::invalidate(std::vector<Node *> *old_nodes, 
                                 bool delete_source = false) {

    std::vector<Node *> to_remove;
    to_remove.push_back(this);

    while (!to_remove.empty()) {
        Node *n = to_remove.back();
        to_remove.pop_back();

        for (auto p = n->parents_.begin(); p != n->parents_.end(); p++) {
            if (p->first->remove_child(n)) {
                to_remove.push_back(p->first);
            }
        }

        if (!n->parents_.empty() || delete_source) {
            old_nodes->push_back(n);
        }

        n->parents_.clear();
        n->children_.clear();

        n->length_ = n->stay_count_ = n->seed_prob_ = 0;
    }
}

void SeedGraph::Node::replace_info(const Node &node) {
    length_ = node.length_;

    event_prob_ = node.event_prob_;
    seed_prob_ = node.seed_prob_;

    stay_count_ = node.stay_count_;
    skip_count_ = node.skip_count_;
    ignore_count_ = node.ignore_count_;

    #ifdef DEBUG_PROB
    dup_node->min_evt_prob_ = node.min_evt_prob_;
    #endif
}

void SeedGraph::Node::update_info() {
    if (parents_.empty()) {
        seed_prob_ = event_prob_;
        stay_count_ = 
        skip_count_ =
        ignore_count_ = 0;

        #ifdef DEBUG_PROB
        min_evt_prob_ = seed_prob_;
        #endif

        return;
    }

    const parent_ptr &parent = parents_.front();
    
    seed_prob_    = parent.first->seed_prob_ + event_prob_;

    stay_count_   = parent.first->stay_count_;
    skip_count_   = parent.first->skip_count_;
    ignore_count_ = parent.first->ignore_count_;

    switch(parent.second) {
        case Type::STAY:
        stay_count_++;
        break;

        case Type::SKIP:
        skip_count_++;
        break;

        case Type::IGNORE:
        ignore_count_++;

        default:
            break;
    }

    #ifdef DEBUG_PROB
    if (parent.first->min_evt_prob_ < event_prob_) {
        min_evt_prob_ = parent.first->min_evt_prob_;
    } else {
        min_evt_prob_ = event_prob_;
    }
    #endif
}


//Returns true if this node should be pruned
//If should_invalidate = false, will always return false
bool SeedGraph::Node::remove_child(Node *child) {

    if (children_.size() == 1) {
        children_.clear();
        return true;
    }

    for (auto c = children_.begin(); c != children_.end(); c++) {
        if (*c == child) {
            children_.erase(c);
            break;
        }
    }

    return false;
}



//SeedGraph::SeedGraph(const KmerModel &model,
//                     const KmerFMI &fmi, 
//                     const NormParams &norm_params, 
//                     int seed_len, int read_len,
//                     double min_event_prob,
//                     double min_seed_prob,
//                     double min_stay_prob,
//                     double max_stay_frac,
//                     const std::string &label)
//    : model_(model),
//      fmi_(fmi),
//      norm_params_(norm_params),
//      seed_length_(seed_len),
//      cur_event_(read_len), 
//      max_stays_( (int) (max_stay_frac * seed_len) ),
//      prev_event_({0, 0, 0}),
//      min_event_prob_(min_event_prob),
//      min_seed_prob_(min_seed_prob),
//      min_stay_prob_(min_stay_prob) {
//}
//

AlnParams::AlnParams(const KmerModel &model,
                     int min_seed_nlen, 
                     int anchor_nlen, 
                     int max_ignores, 
                     int max_skips,
                     double max_stay_frac,
                     double min_anchor_evpr,
                     double min_extend_evpr,
                     double min_seed_pr,
                     double min_stay_pr) 
        : model_(model),
          max_ignores_(max_ignores),
          max_skips_(max_skips),
          max_stay_frac_(max_stay_frac),
          min_anchor_evpr_(min_anchor_evpr),
          min_extend_evpr_(min_extend_evpr),
          min_seed_pr_(min_seed_pr),
          min_stay_pr_(min_stay_pr) {

    anchor_rlen_ = nucl_to_events(anchor_nlen);
    graph_elen_ = get_graph_len(min_seed_nlen);

    std::cerr << "Graph len: " << graph_elen_ << "\n";

}

int AlnParams::nucl_to_events(int n) {
    return n - model_.kmer_len() + 1;
}

int AlnParams::get_graph_len(int seed_nlen) {
    return (nucl_to_events(seed_nlen) / (1.0 - max_stay_frac_)) + max_ignores_;
}

SeedGraph::SeedGraph(const KmerFMI &fmi, 
                     const AlnParams &ap,
                     const std::string &label)
    : fmi_(fmi),
      params_(ap),
      label_(label) {
    timer.reset();
}

SeedGraph::~SeedGraph() {
    reset();
}

void SeedGraph::new_read(int read_len, const NormParams &params) {
    reset();
    norm_params_ = params;
    cur_event_ = read_len;
    prev_event_ = {0, 0, 0};
}

void SeedGraph::reset() {
    for (auto p = prev_nodes_.begin(); p != prev_nodes_.end(); p++) {
        if (!p->first->parents_.empty())
            p->first->invalidate(&old_nodes_);
    }

    prev_nodes_.clear();

    for (auto ss = sources_.begin(); ss != sources_.end(); ss++) {
        for (auto s = ss->begin(); s != ss->end(); s++) {
            delete *s;
        }
    }

    sources_.clear();

    for (int i = 0; i < old_nodes_.size(); i++) {
        delete old_nodes_[i];
    }

    old_nodes_.clear();

    for (auto probs = event_kmer_probs_.begin(); 
         probs != event_kmer_probs_.end(); probs++) {
        delete [] *probs;
    }
    timer.reset();

    event_kmer_probs_.clear();
}

std::vector<Result> SeedGraph::add_event(Event e, std::ostream &out) {
    //Update event index
    cur_event_--;

    if (!params_.model_.event_valid(e)) {
        std::cerr << "Error: event " << cur_event_
                  << " invalid - this might cause some problems\n";
    }

    //Compare this event with previous to see if it should be a stay

    double cur_stay_prob = params_.model_.get_stay_prob(e, prev_event_);

    bool stay = true;//cur_stay_prob >= min_stay_prob_;


    Timer t;

    //Calculate and store kmer match probs for this event
    event_kmer_probs_.push_front(new double[params_.model_.kmer_count()]);
    double *kmer_probs = event_kmer_probs_.front();
    for (int kmer = 0; kmer < params_.model_.kmer_count(); kmer++)
        kmer_probs[kmer] = params_.model_.event_match_prob(e, kmer, norm_params_);

    //DEBUG(cur_event_ << "\t" << t.lap() << "\t");
    
    //Where this event's sources will be stored
    sources_.push_front(std::list<Node *>());
    
    bool event_skipped = false;

    double prob;

    //Find neighbors of previous nodes
    for (auto p = prev_nodes_.begin(); p != prev_nodes_.end(); p++) {
        Range prev_range = p->second;
        Node *prev_node = p->first;

        //Get probability for stay neighbor
        int prev_kmer = prev_node->kmer_;
        prob = kmer_probs[prev_kmer];
        
        int neighbor_count = 0;

        double evpr_thresh;

        if (prev_node->match_len() < params_.anchor_rlen_ &&
            prev_range.length() > 1) { 
            evpr_thresh = params_.min_anchor_evpr_;
        } else {
            evpr_thresh = params_.min_extend_evpr_;
        }

        if (stay && prob >= evpr_thresh) {
            neighbor_count += 1;

            Node next_node(prev_node, prev_kmer, prob, Node::Type::STAY);

            add_child(prev_range, next_node);
        }
        
        //Find next possible kmers
        auto neighbor_itr = params_.model_.get_neighbors(prev_kmer);
        std::list<mer_id> next_kmers;

        for (auto n = neighbor_itr.first; 
             n != neighbor_itr.second; n++) {


            if(kmer_probs[*n] >= evpr_thresh) {
                next_kmers.push_back(*n);
            } 
        }

        //Find ranges FM index for those next kmers
        std::list<Range> next_ranges = fmi_.get_neigbhors(prev_range, next_kmers);

        auto next_kmer = next_kmers.begin(); 
        auto next_range = next_ranges.begin();
        

        //Add all the neighbors that were found
        while (next_kmer != next_kmers.end()) {
            prob = kmer_probs[*next_kmer];

            neighbor_count += next_range->is_valid();

            //std::cout << next_range->start_ << " " << next_range->end_ << "\n";
            Node next_node(prev_node, *next_kmer, prob, Node::Type::MATCH);
            add_child(*next_range, next_node);

            next_kmer++; 
            next_range++;
        }

        if (neighbor_count < 2 && //Probably ok? Maybe add param?

            //Maybe mark if previously a full seed, so only happens when extending
            prev_node->seed_len() == params_.graph_elen_ - 1 && 
            prev_range.length() == 1 && //Yes
            prev_node->ignore_count_ <= params_.max_ignores_) { //max_ignores (frac?)

            prob = prev_node->seed_prob_ / prev_node->length_; //could be better
            
            Node next_node(prev_node, prev_kmer, prob, Node::Type::IGNORE);
            add_child(prev_range, next_node);
        }
    }

    //DEBUG(t.lap() << "\t");

    //Find sources
    for (mer_id kmer = 0; kmer < params_.model_.kmer_count(); kmer++) {
        prob = kmer_probs[kmer];
        if (prob >= params_.min_anchor_evpr_) {
            Range next_range = fmi_.get_full_range(kmer);

            if (next_range.is_valid()) {
                Node next_node(kmer, prob);

                //Will split sources if they intersect existing nodes
                add_sources(next_range, next_node);             
            }
        }
    }

    //DEBUG(t.lap() << "\t");

    //Clear prev_nodes, pruning any leaves
    while (!prev_nodes_.empty()) {
        auto p = prev_nodes_.begin();
        Node *prev_node = p->first;

        if (prev_node->children_.empty()) {
            prev_node->invalidate(&old_nodes_);
        }

        prev_nodes_.erase(p);
    }
 
    //Move next_nodes to prev
    while (!next_nodes_.empty()) {
        auto n = next_nodes_.begin();
        prev_nodes_[n->second] = n->first;
        next_nodes_.erase(n);
    }

    //DEBUG(t.lap() << "\t");

    auto r = pop_seeds(out);

    prev_event_ = e;

    return r;
}

int SeedGraph::add_sources(const Range &range, const Node &node) {

    //Find closest existing node
    auto lb = next_nodes_.lower_bound(range);
    
    //Find range of existing nodes intersecting the new node
    auto start = lb;
    while (start != next_nodes_.begin() 
           && std::prev(start)->first.intersects(range)){
        start--;
    }

    auto end = lb;
    while (end != next_nodes_.end() && end->first.intersects(range)) {
        end++;
    }

    //No nodes intersect new node, just add it
    if (start == next_nodes_.end() || (start == lb && end == lb)) {

        if (old_nodes_.empty()) {
            next_nodes_[range] = new Node(node);
        } else {
            *(old_nodes_.back()) = node;
            next_nodes_[range] = old_nodes_.back();
            old_nodes_.pop_back();
        }

        sources_.front().push_back(next_nodes_[range]);
        return 1;
    }
    

    //Split range to not include intersecting ranges
    std::list<Range> split_ranges;

    Range rr(range); //Copy full range

    for (auto pr = start; pr != end; pr++) {

        //After this, rr stores right segment, lr stores left
        Range lr = rr.split_range(pr->first); //Left range

        //Add left range if it exists
        if (lr.is_valid()) {
            split_ranges.push_back(lr);
        }
        
        //If right range invalid, no more splitting can be done
        if (!rr.is_valid())
            break;
    }

    //Add last range if there is one
    if (rr.is_valid()) 
        split_ranges.push_back(rr);

    //Add a new node for every split range
    for (auto r = split_ranges.begin(); r != split_ranges.end(); r++) {
        if (old_nodes_.empty()) {
            next_nodes_[*r] = new Node(node);
        } else {
            *(old_nodes_.back()) = node;
            next_nodes_[*r] = old_nodes_.back();
            old_nodes_.pop_back();
        }

        sources_.front().push_back(next_nodes_[*r]);
    }

    return split_ranges.size();
}

SeedGraph::Node *SeedGraph::add_child(Range &range, Node &node) {
    
    if (!node.is_valid() || !range.is_valid()) {
        return NULL;
    }
    
    //Parent of the node being added
    parent_ptr new_parent = node.parents_.front();

    //Find closest node >= the node being added
    auto lb = next_nodes_.lower_bound(range);

    //Node range hasn't been added yet
    if (lb == next_nodes_.end() || !lb->first.same_range(range)) {

        //Allocate memory for a new node
        Node *new_node;
        if (old_nodes_.empty()) {
            new_node = new Node(node);
        } else {
            *(old_nodes_.back()) = node;
            new_node = old_nodes_.back();
            old_nodes_.pop_back();
        }

        //Store it with it's range
        next_nodes_.insert(lb, std::pair<Range, Node *>(range, new_node));

        //Update the parent
        new_parent.first->children_.push_front(new_node);

        //Created one node
        return new_node;
    }
    
    //Node associated with same range
    Node *dup_node = lb->second;
    auto dup_parent = dup_node->parents_.begin(); 

    //This node is for a longer seed
    if (dup_node->seed_len() < node.seed_len()) {

        //Copy seed info
        dup_node->replace_info(node);

        //Update parents and children
        dup_node->parents_.push_front(new_parent);
        new_parent.first->children_.push_front(dup_node);

        //No new nodes created
        return dup_node;

    } else if (dup_node->seed_len() == node.seed_len()) {

        if(node.better_than(dup_node)) {
        //if ( dup_node->ignore_count_ > node.ignore_count_ 
        //    || (dup_node->ignore_count_ == node.ignore_count_ && 
        //        dup_node->seed_prob_ < node.seed_prob_)) {

            //If dup_parent has no children in the end, will be invalidated (in add_event)
            dup_parent->first->remove_child(dup_node);

            //Erase old parent 
            dup_node->parents_.erase(dup_parent);

            //Copy seed info
            dup_node->replace_info(node);

            //Insert new parent in old parent's place
            dup_node->parents_.push_front(new_parent);
            new_parent.first->children_.push_front(dup_node);
        }

        return dup_node;
    }
   
    //Find where to place the new parent pointer 
    while (dup_parent != dup_node->parents_.end()) {

        //dup_parent is a longer seed
        if (new_parent.first->seed_len() < dup_parent->first->seed_len()) {
            dup_parent++;
            continue; //not breaking yet

        //dup_parent is the longest branch shorter than new_parent
        } else if (new_parent.first->seed_len() > dup_parent->first->seed_len())  {
            
            //Add parent in the correct location
            dup_node->parents_.insert(dup_parent, new_parent);
            new_parent.first->children_.push_front(dup_node);

        //dup_parent is same length as new_parent

        //and new node has higher probability
        } else if (new_parent.first->better_than(dup_parent->first)) {

            //If dup_parent has no children in the end, will be invalidated (in add_event)
            dup_parent->first->remove_child(dup_node);

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
        dup_node->parents_.push_back(new_parent);
        new_parent.first->children_.push_front(dup_node);
    }

    return dup_node;
}



std::vector<Result> SeedGraph::pop_seeds(std::ostream &out) { 

    std::vector<Result> results;

    //Is there a long enough seed?
    if (sources_.size() < params_.graph_elen_)
        return results;

    std::list<Node *> &aln_ends = sources_.back(),
                      &next_sources = *std::next(sources_.rbegin());

    while (!aln_ends.empty()) {
        Node *aln_en = aln_ends.front();

        if (!aln_en->is_valid()) {
            aln_ends.pop_front();
            aln_en->invalidate(&old_nodes_, true);
            //delete aln_en;
            continue;
        }

        std::list<Node *> to_visit;
        for (auto c = aln_en->children_.begin(); c != aln_en->children_.end(); c++) {
            (*c)->parents_.clear();
            next_sources.push_back(*c);
            to_visit.push_back(*c);
        }
        
        double prev_len = aln_en->seed_len();
        auto kmer_probs = event_kmer_probs_.rbegin();

        while (!to_visit.empty()) {


            Node *n = to_visit.front();
            to_visit.pop_front();

            if (n->seed_len() != prev_len) {
                prev_len = n->seed_len();
                kmer_probs++;
            }

            //Check if parents are compatible
            if (n->parents_.size() > 1) {
                auto seed_parent = n->parents_.begin(), //where we came from
                     next_parent = std::next(seed_parent); //seed with closest length

                //Parents now have the same length
                if (seed_parent->first->seed_len() == next_parent->first->seed_len()) {
                    Node *to_erase;

                    //Erase where we came from
                    if (next_parent->first->better_than(seed_parent->first)) {
                    //if (seed_parent->first->seed_prob_ < next_parent->first->seed_prob_) {
                        to_erase = seed_parent->first; 
                        n->parents_.erase(seed_parent);

                    //Erase other seed branch
                    } else {
                        to_erase = next_parent->first;
                        n->parents_.erase(next_parent);
                    }

                    //Remove old parent
                    if (to_erase->remove_child(n)) {
                        to_erase->invalidate(&old_nodes_);
                    }
                }

            }

            if (n->seed_len() == params_.graph_elen_) {

                if (n->should_report(params_)) {
                    Range range = prev_nodes_[n];
                    Result r(cur_event_, params_.graph_elen_, n->seed_prob_ / n->seed_len());

                    #ifdef DEBUG_PROB
                    r.min_evt_prob_ = n->min_evt_prob_ ;
                    #endif

                    for (int s = range.start_; s <= range.end_; s++) {
                        r.set_ref_range(fmi_.suffix_ar_[s], n->match_len());
                        results.push_back(r);

                        out << label_ << "\t" << timer.get() << "\t";
                        r.print(out);
                    }
                }

            } else {
                for (auto c = n->children_.begin(); c != n->children_.end(); c++) {
                    to_visit.push_back(*c);
                }
            }

            n->length_--;
            
            n->update_info();

        }
        
        aln_ends.pop_front();
        aln_en->invalidate(&old_nodes_, true);
        //delete aln_en;
    }

    delete [] event_kmer_probs_.back();
    event_kmer_probs_.pop_back();
    sources_.pop_back();

    return results;

}

void SeedGraph::print_graph(bool verbose) {
    std::set< std::pair<int, Node *> > to_visit;
    auto event_sources = sources_.rbegin();

    std::cout << "== printing graph at event " << cur_event_ << " ==\n"; 

    for (auto s = event_sources->begin(); s != event_sources->end(); s++) 
        to_visit.insert( std::pair<int, Node *> (0, *s) );

    bool iter_sources = true;

    int source_count = 0;

    std::vector<int> source_counts(sources_.size(), 0),
                     linear_counts(sources_.size(), 0),
                     branch_counts(sources_.size(), 0),
                     invalid_counts(sources_.size(), 0);

    while (!to_visit.empty()) {

        if (iter_sources) {
            source_count++;
            event_sources++;
        }

        //TODO: move this outside
        if (event_sources != sources_.rend()) {
            for (auto s = event_sources->begin(); s != event_sources->end(); s++) 
                to_visit.insert(std::pair<int, Node *>(source_count, *s));
        } else {
            iter_sources = false;
        }

        auto p = to_visit.begin();
        int event = p->first;
        Node *n = p->second;
        to_visit.erase(p);

        if (!n->is_valid()) {
            invalid_counts[event]++;
        } else if (n->parents_.empty()) {
            source_counts[event]++;
        } else if (n->children_.size() > 1 || n->parents_.size() > 1) {
            branch_counts[event]++;
        } else {
            linear_counts[event]++;
        }
        
        if (verbose) {
            std::cout << n << " " << event << " " 
                      << n->seed_len() << " " 
                      << n->parents_.size();
        }

        for (auto c = n->children_.begin(); c != n->children_.end(); c++) {
            if (verbose) {
                std::cout << " " << *c;
            }
            to_visit.insert(std::pair<int, Node *>(event + 1, *c));
        }

        if (verbose) {
            std::cout << "\n";
        }
    }

    int source_total = 0;
    std::cout << "== source counts:";
    for (unsigned int i = 0; i < source_counts.size(); i++) {
        std::cout << "\t" << source_counts[i];
        source_total += source_counts[i];
    }
    std::cout << "\t (" << source_total << ")\t==\n";

    int branch_total = 0;
    std::cout << "== branch counts:";
    for (unsigned int i = 0; i < branch_counts.size(); i++) {
        std::cout << "\t" << branch_counts[i];
        branch_total += branch_counts[i];
    }
    std::cout << "\t (" << branch_total << ")\t==\n";


    int linear_total = 0;
    std::cout << "== linear counts:";
    for (unsigned int i = 0; i < linear_counts.size(); i++) {
        std::cout << "\t" << linear_counts[i];
        linear_total += linear_counts[i];
    }
    std::cout << "\t (" << linear_total << ")\t==\n";

    int invalid_total = 0;
    std::cout << "== invalid counts:";
    for (unsigned int i = 0; i < invalid_counts.size(); i++) {
        std::cout << "\t" << invalid_counts[i];
        invalid_total += invalid_counts[i];
    }
    std::cout << "\t (" << invalid_total << "\t==\n";

    std::cout << "== full total: " 
              << (source_total + branch_total + linear_total + invalid_total) 
              << "\t==\n";

    std::cout << "== end ==\n";
}

Result::Result(int read_start, int seed_len, double prob, int ref_start, int ref_end) 
    : read_range_( Range(read_start, read_start + seed_len - 1) ),
      ref_range_(Range(ref_start, ref_end)),
      seed_prob_(prob) {}

void Result::set_ref_range(int start, int length) {
    ref_range_.start_ = start;
    ref_range_.end_ = start + length - 1;
}


void Result::print(std::ostream &out) {
    out << read_range_.start_ << "-" << read_range_.end_ << "\t"
              << ref_range_.start_ << "-" << ref_range_.end_ << "\t"
              << seed_prob_;
              

    #ifdef DEBUG_PROB
    out << " " << min_evt_prob_;
    #endif

    out << "\n";
}
