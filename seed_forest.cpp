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
#include <sdsl/suffix_arrays.hpp>
#include "timer.hpp"
#include "seed_forest.hpp"
#include "fmi.hpp"

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)
//#define DEBUG_PATHS
//#define DEBUG_RANGES


//Source constructor
SeedForest::Node::Node(mer_id kmer, double prob)
    : length_(1),
      stay_count_(0),
      skip_count_(0),
      ignore_count_(0),
      consec_stays_(0),
      kmer_(kmer),
      event_prob_(prob),
      seed_prob_(prob)

      #ifdef DEBUG_NODES
      , full_length_(1)
      #endif

      #ifdef DEBUG_PROB
      , min_evt_prob_(prob)
      #endif
      
       {}

//Child constructor
SeedForest::Node::Node(Node *parent, mer_id kmer, 
                      double prob, Type type) 
    :
      kmer_(kmer),
      event_prob_(prob) {

    #ifdef DEBUG_NODES
    full_length_ = parent->full_length_ + 1;
    #endif
    
    parents_.push_front(parent_ptr(parent, type));
    update_info();
}

//Creates invalid node
SeedForest::Node::Node()
    : length_(0),
      stay_count_(0),
      skip_count_(0),
      ignore_count_(0),
      consec_stays_(0),
      kmer_(0),
      event_prob_(0),
      seed_prob_(0) {}

//Copy constructor
SeedForest::Node::Node(const Node &s)
    : length_(s.length_),
      stay_count_(s.stay_count_),
      skip_count_(s.skip_count_),
      ignore_count_(s.ignore_count_),
      consec_stays_(s.consec_stays_),
      kmer_(s.kmer_),
      event_prob_(s.event_prob_),
      seed_prob_(s.seed_prob_),

      #ifdef DEBUG_PROB
      min_evt_prob_(s.min_evt_prob_),
      #endif

      #ifdef DEBUG_NODES
      full_length_(s.full_length_),
      #endif

      parents_(s.parents_),
      children_(s.children_) {}
    
bool SeedForest::Node::is_valid() {
    return length_ > 0;
}

size_t SeedForest::Node::seed_len() {
    return length_;
}

size_t SeedForest::Node::match_len() {
    return seed_len() - stay_count_ - ignore_count_ + skip_count_;
}

double SeedForest::Node::mean_prob() const {
    return seed_prob_ / length_;
}

//TODO: Really not a very good name
bool SeedForest::Node::better_than(const Node *node) {
    return node->ignore_count_ > ignore_count_ ||
            (node->ignore_count_ == ignore_count_ && 
             node->mean_prob() < mean_prob());
}

bool SeedForest::Node::should_report(const AlnParams &params) {
    return parents_.front().second == Node::Type::MATCH && 
           stay_count_ <= params.max_stay_frac_ * length_ && 
           seed_prob_ >= params.min_seed_pr_ * length_;
}

void SeedForest::Node::invalidate(std::vector<Node *> *old_nodes, 
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

void SeedForest::Node::replace_info(const Node &node) {
    length_ = node.length_;

    event_prob_ = node.event_prob_;
    seed_prob_ = node.seed_prob_;

    stay_count_ = node.stay_count_;
    skip_count_ = node.skip_count_;
    ignore_count_ = node.ignore_count_;
    consec_stays_ = node.consec_stays_;

    #ifdef DEBUG_PROB
    dup_node->min_evt_prob_ = node.min_evt_prob_;
    #endif

    #ifdef DEBUG_NODES
    full_length_ = node.full_length_;
    #endif
}

void SeedForest::Node::update_info() {
    if (parents_.empty()) {
        seed_prob_ = event_prob_;
        stay_count_ = 
        skip_count_ =
        ignore_count_ = 
        consec_stays_ = 0;

        length_ = 1;

        #ifdef DEBUG_PROB
        min_evt_prob_ = seed_prob_;
        #endif

        return;
    }

    const parent_ptr &parent = parents_.front();
    
    seed_prob_    = parent.first->seed_prob_ + event_prob_;
    
    length_       = parent.first->length_ + 1;
    stay_count_   = parent.first->stay_count_;
    skip_count_   = parent.first->skip_count_;
    ignore_count_ = parent.first->ignore_count_;
    consec_stays_ = 0;

    switch(parent.second) {
        case Type::STAY:
        stay_count_++;
        consec_stays_ = parent.first->consec_stays_ + 1;
        break;

        case Type::SKIP:
        skip_count_++;
        break;

        case Type::IGNORE:
        ignore_count_++;

        default:
            break;
    }

    #ifdef DEBUG_NODES
    full_length_ = parent.first->full_length_ + 1;
    #endif

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
bool SeedForest::Node::remove_child(Node *child) {

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

void SeedForest::Node::print() const {
    //std::cout << length_ << "\n";
              //<< full_length_ << "\n";
              //<< stay_count_ << "\t"
              ////<< skip_count_ << "\t"
              ////<< ignore_count_ << "\t"
              //<< (parents_.empty() ? Type::MATCH : parents_.front().second) << "\t"
              //<< parents_.size() << "\t"
              //<< children_.size() << "\t"
              //<< kmer_ << "\t"
              //<< event_prob_ << "\t"
              //<< seed_prob_ << "\n";
}


AlnParams::AlnParams(const KmerModel &model,
                     unsigned int min_seed_nlen, 
                     unsigned int anchor_nlen, 
                     unsigned int max_ignores, 
                     unsigned int max_skips,
                     unsigned int max_consec_stay,
                     double max_stay_frac,
                     double min_anchor_evpr,
                     //double min_extend_evpr,
                     std::vector<unsigned int>    expr_lengths,
                     std::vector<double> expr_probs,
                     double min_seed_pr,
                     double min_stay_pr)
        : model_(model),
          max_ignores_(max_ignores),
          max_skips_(max_skips),
          max_consec_stay_(max_consec_stay),
          max_stay_frac_(max_stay_frac),
          min_anchor_evpr_(min_anchor_evpr),
          min_seed_pr_(min_seed_pr),
          min_stay_pr_(min_stay_pr),
          //min_extend_evpr_(min_extend_evpr),
          expr_lengths_(expr_lengths),
          expr_probs_(expr_probs) {

    anchor_rlen_ = nucl_to_events(anchor_nlen);
    graph_elen_ = get_graph_len(min_seed_nlen);

    std::cerr << "Graph len: " << graph_elen_ << "\n";

}

unsigned int AlnParams::nucl_to_events(unsigned int n) {
    return n - model_.kmer_len() + 1;
}

unsigned int AlnParams::get_graph_len(unsigned int seed_nlen) {
    return (nucl_to_events(seed_nlen) / (1.0 - max_stay_frac_)) + max_ignores_;
}

SeedForest::SeedForest(const FMI &fmi, 
                     const AlnParams &ap,
                     const std::string &label)
    : fmi_(fmi),
      params_(ap),
      label_(label) {
    timer.reset();

    kmer_ranges_ = new Range[params_.model_.kmer_count()];
    for (mer_id k = 0; k < params_.model_.kmer_count(); k++) {
        Range r = fmi_.get_full_range(params_.model_.get_last_base(k));
        for (size_t i = params_.model_.kmer_len()-2; 
             i < params_.model_.kmer_len(); i--) {
            r = fmi_.get_neighbor(r, params_.model_.get_base(k, i));
        }
        kmer_ranges_[k] = r;
    }
}

SeedForest::~SeedForest() {
    delete[] kmer_ranges_;
    reset();
}

void SeedForest::new_read(size_t read_len) {
    reset();
    cur_event_ = read_len;
    prev_event_ = {0, 0, 0};
}

void SeedForest::reset() {
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

    for (size_t i = 0; i < old_nodes_.size(); i++) {
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

std::vector<Result> SeedForest::add_event(Event e, std::ostream &out) {
    //Update event index
    cur_event_--;

    Timer t;

    //Calculate and store kmer match probs for this event
    event_kmer_probs_.push_front(new double[params_.model_.kmer_count()]);
    double *kmer_probs = event_kmer_probs_.front();
    for (unsigned int kmer = 0; kmer < params_.model_.kmer_count(); kmer++) {
        kmer_probs[kmer] = params_.model_.event_match_prob(e, kmer);
    }

    //Where this event's sources will be stored
    sources_.push_front(std::list<Node *>());
    
    double prob;

    t.reset();

    //Find neighbors of previous nodes
    for (auto p = prev_nodes_.begin(); p != prev_nodes_.end(); p++) {
        Range prev_range = p->second;
        Node *prev_node = p->first;

        #ifdef DEBUG_RANGES
        //if (prev_range.length() <= 200) {
        //}
        #endif

        //Get probability for stay neighbor
        mer_id prev_kmer = prev_node->kmer_;
        prob = kmer_probs[prev_kmer];
        //rank = kmer_ranks[prev_kmer];
        
        size_t neighbor_count = 0;

        double evpr_thresh = 0;
        for (size_t i = 0; i < params_.expr_lengths_.size(); i++) {
            if (prev_range.length() <= params_.expr_lengths_[i]) {
                evpr_thresh = params_.expr_probs_[i];
                break;
            }
        }

        if (evpr_thresh == 0) {
            evpr_thresh = params_.min_anchor_evpr_;
        }

        if (prev_node->consec_stays_ < params_.max_consec_stay_ && prob >= evpr_thresh) { //&& rank <= evpr_thresh ) { 
            neighbor_count += 1;

            Node next_node(prev_node, prev_kmer, prob, Node::Type::STAY);

            #ifdef DEBUG_RANGES
            //if (prev_range.length() <= 200) {
                std::cout << "\t" << prev_range.length();
            //}
            #endif

            add_child(prev_range, next_node);
        }
        
        //Find next possible kmers
        auto neighbor_itr = params_.model_.get_neighbors(prev_kmer);
        std::list<mer_id> next_kmers;

        for (auto n = neighbor_itr.first; n != neighbor_itr.second; n++) {
            if(kmer_probs[*n] >= evpr_thresh) {
            //if(kmer_ranks[*n] <= evpr_thresh) {
                next_kmers.push_back(*n);
            } 
        }

        //Find ranges FM index for those next kmers
        //std::list<Range> next_ranges = fmi_.get_neigbhors(prev_range, next_kmers);

        auto next_kmer = next_kmers.begin(); 
        //auto next_range = next_ranges.begin();
        

        //Add all the neighbors that were found
        while (next_kmer != next_kmers.end()) {
            prob = kmer_probs[*next_kmer];
            //rank = kmer_ranks[*next_kmer];

            base_t next_base = params_.model_.get_first_base(*next_kmer);
            Range next_range = fmi_.get_neighbor(prev_range, next_base);

            neighbor_count += next_range.is_valid();

            #ifdef DEBUG_RANGES
            //if (prev_range.length() <= 200 && next_range.is_valid()) {
            if (next_range.is_valid()) {
                std::cout << "\t" << next_range.length();

            }
            #endif

            Node next_node(prev_node, *next_kmer, prob, Node::Type::MATCH);
            add_child(next_range, next_node);

            next_kmer++; 
            //next_range++;
        }

        if (neighbor_count < 2 && //Probably ok? Maybe add param?

            //Maybe mark if previously a full seed, so only happens when extending
            prev_node->seed_len() == params_.graph_elen_ - 1 && 
            prev_range.length() == 1 && //Yes
            prev_node->ignore_count_ < params_.max_ignores_) { //max_ignores (frac?)

            prob = prev_node->seed_prob_ / prev_node->length_; //could be better
            
            Node next_node(prev_node, prev_kmer, prob, Node::Type::IGNORE);
            add_child(prev_range, next_node);
        }
        #ifdef DEBUG_RANGES
        //if (prev_range.length() <= 200) {
            std::cout << "\n";
        //}
        #endif
    }

    //if (cur_event_ < 5000) {
    //    std::cout << timer.lap() << std::endl;
    //}

    //DEBUG(t.lap() << "\t");

    //Find sources
    for (mer_id kmer = 0; kmer < params_.model_.kmer_count(); kmer++) {
        prob = kmer_probs[kmer];
        //rank = kmer_ranks[kmer];
        if (prob >= params_.min_anchor_evpr_) {
        //if (rank <= params_.min_anchor_evpr_) {
            //Range next_range = fmi_.get_kmer_range(params_.model_.id_to_kmer(kmer));
            Range next_range = kmer_ranges_[kmer];

            if (next_range.is_valid()) {
                Node next_node(kmer, prob);
                //std::cout << kmer << "\t" << prob << "\n";

                //Will split sources if they intersect existing nodes
                add_sources(next_range, next_node);             
            }
        }
    }

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

    //if (cur_event_ < 5000) {
    //    std::cout << timer.lap() << std::endl;
    //}

    auto r = pop_seeds(out);

    //if (cur_event_ < 5000) {
    //    std::cout << timer.lap() << std::endl;
    //}

    prev_event_ = e;

    return r;
}

size_t SeedForest::add_sources(const Range &range, const Node &node) {

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


        //#ifdef DEBUG_NODES
        //std::cout << range.start_ << "\t" << range.end_ << "\t";
        //node.print();
        //#endif

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
        //#ifdef DEBUG_NODES
        //std::cout << r->start_ << "\t" << r->end_ << "\t";
        //node.print();
        //#endif

        sources_.front().push_back(next_nodes_[*r]);
    }

    return split_ranges.size();
}

SeedForest::Node *SeedForest::add_child(Range &range, Node &node) {
    
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

    if(node.better_than(dup_node)) {
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



std::vector<Result> SeedForest::pop_seeds(std::ostream &out) { 

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

        #ifdef DEBUG_PATHS
        std::list< std::vector<Node::Type> > paths;
        #endif 

        std::list<Node *> to_visit;

        for (auto c = aln_en->children_.begin(); c != aln_en->children_.end(); c++) {

            #ifdef DEBUG_PATHS
            std::vector<Node::Type> p(params_.graph_elen_);
            p[0] = Node::Type::MATCH;
            p[1] = (*c)->parents_.front().second;
            paths.push_back(p);
            #endif 

            (*c)->parents_.clear();
            next_sources.push_back(*c);
            to_visit.push_back(*c);

        }
        
        double prev_len = aln_en->seed_len();
        auto kmer_probs = event_kmer_probs_.rbegin();

        while (!to_visit.empty()) {
    
            Node *n = to_visit.front();
            to_visit.pop_front();

            #ifdef DEBUG_PATHS
            std::vector<Node::Type> p = paths.front();
            paths.pop_front();

            if (n->seed_len() > 2) {
                p[n->seed_len()-1] = n->parents_.front().second;
            }
            #endif 

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

                    n->update_info();
                }

            }

            if (n->seed_len() == params_.graph_elen_) {

                Range range = prev_nodes_[n];

                //TODO: max repeat parameter
                if (n->should_report(params_) && range.length() < 100) { 
                    Result r(cur_event_, params_.graph_elen_, n->seed_prob_ / n->seed_len());

                    #ifdef DEBUG_PROB
                    r.min_evt_prob_ = n->min_evt_prob_ ;
                    #endif

                    for (unsigned int s = range.start_; s <= range.end_; s++) {
                        r.set_ref_range(fmi_.sa(s), n->match_len());
                        results.push_back(r);

                        out << label_ << "\t" << timer.get() << "\t";

                        #ifdef DEBUG_PATHS
                        for (size_t i = p.size()-1; i < p.size(); i--) {
                            switch (p[i]) {
                                case Node::Type::MATCH:
                                out << "M";
                                break;
                                case Node::Type::STAY:
                                out << "S";
                                break;
                                case Node::Type::SKIP:
                                out << "K";
                                break;
                                case Node::Type::IGNORE:
                                out << "I";
                                break;
                                default:
                                out << "?";
                                break;
                            }
                        }
                        out << "\t";
                        #endif
                        r.print(out);
                    }
                }

            } else {

                #ifdef DEBUG_PATHS
                if (n->children_.size() > 0) { //Should always be true?
                    paths.push_back(p);

                    for (unsigned int i = 1; i < n->children_.size(); i++) {
                        paths.push_back(std::vector<Node::Type>(p));
                    }   
                }
                #endif

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

void SeedForest::print_graph(bool verbose) {
    std::set< std::pair<unsigned int, Node *> > to_visit;
    auto event_sources = sources_.rbegin();

    std::cout << "== printing graph at event " << cur_event_ << " ==\n"; 

    for (auto s = event_sources->begin(); s != event_sources->end(); s++) 
        to_visit.insert( std::pair<unsigned int, Node *> (0, *s) );

    bool iter_sources = true;

    unsigned int source_count = 0;

    std::vector<unsigned int> source_counts(sources_.size(), 0),
                     linear_counts(sources_.size(), 0),
                     branch_counts(sources_.size(), 0),
                     multiparent_counts(sources_.size(), 0),
                     invalid_counts(sources_.size(), 0);

    while (!to_visit.empty()) {

        if (iter_sources) {
            source_count++;
            event_sources++;
        }

        //TODO: move this outside
        if (event_sources != sources_.rend()) {
            for (auto s = event_sources->begin(); s != event_sources->end(); s++) 
                to_visit.insert(std::pair<unsigned int, Node *>(source_count, *s));
        } else {
            iter_sources = false;
        }

        auto p = to_visit.begin();
        unsigned int event = p->first;
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

        if (n->parents_.size() > 1) {
            multiparent_counts[event]++;
        }
        
        if (verbose) {
            std::cout << n << " " << event << " " 
                      << n->seed_len() << " " 
                      << n->parents_.size() << " ";
            
            if (n->parents_.empty()) {
                std::cout << -1;
            } else {
                std::cout << n->parents_.front().second;
            }
        }

        for (auto c = n->children_.begin(); c != n->children_.end(); c++) {
            if (verbose) {
                std::cout << " " << *c;
            }
            to_visit.insert(std::pair<unsigned int, Node *>(event + 1, *c));
        }

        if (verbose) {
            std::cout << "\n";
        }
    }

    unsigned int source_total = 0;
    std::cout << "== source counts:";
    for (unsigned int i = 0; i < source_counts.size(); i++) {
        std::cout << "\t" << source_counts[i];
        source_total += source_counts[i];
    }
    std::cout << "\t (" << source_total << ")\t==\n";

    unsigned int branch_total = 0;
    std::cout << "== branch counts:";
    for (unsigned int i = 0; i < branch_counts.size(); i++) {
        std::cout << "\t" << branch_counts[i];
        branch_total += branch_counts[i];
    }
    std::cout << "\t (" << branch_total << ")\t==\n";


    unsigned int linear_total = 0;
    std::cout << "== linear counts:";
    for (unsigned int i = 0; i < linear_counts.size(); i++) {
        std::cout << "\t" << linear_counts[i];
        linear_total += linear_counts[i];
    }
    std::cout << "\t (" << linear_total << ")\t==\n";

    unsigned int invalid_total = 0;
    std::cout << "== invalid counts:";
    for (unsigned int i = 0; i < invalid_counts.size(); i++) {
        std::cout << "\t" << invalid_counts[i];
        invalid_total += invalid_counts[i];
    }
    std::cout << "\t (" << invalid_total << ")\t==\n";

    unsigned int multiparent_total = 0;
    std::cout << "== multiparent counts:";
    for (unsigned int i = 0; i < multiparent_counts.size(); i++) {
        std::cout << "\t" << multiparent_counts[i];
        multiparent_total += multiparent_counts[i];
    }
    std::cout << "\t (" << multiparent_total << ")\t==\n";

    std::cout << "== total nodes: " 
              << (source_total + branch_total + linear_total + invalid_total) 
              << "\t==\n";

    std::cout << "== end ==\n";
}

Result::Result(unsigned int read_start, 
              unsigned  int seed_len, 
              double prob, 
              unsigned int ref_start, 
              unsigned int ref_end) 
    : read_range_( Range(read_start, read_start + seed_len - 1) ),
      ref_range_(Range(ref_start, ref_end)),
      seed_prob_(prob) {}

void Result::set_ref_range(unsigned int start, unsigned int length) {
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
