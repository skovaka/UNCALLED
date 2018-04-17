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
#include "forest_aligner.hpp"
#include "fmi.hpp"

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)
//#define DEBUG_PATHS
//#define DEBUG_RANGES

#define MAX_CHILDREN 5

//Source constructor
ForestAligner::Node::Node(Kmer kmer, double prob)
    : length_(1),
      stay_count_(0),
      skip_count_(0),
      ignore_count_(0),
      consec_stays_(0),
      child_count_(0),
      kmer_(kmer),
      event_prob_(prob),
      seed_prob_(prob),
      type_(Type::MATCH),
      parent_(NULL),
      children_(NULL)

      #ifdef DEBUG_NODES
      , full_length_(1)
      #endif

      #ifdef DEBUG_PROB
      , min_evt_prob_(prob)
      #endif
      
       {}

//Child constructor
ForestAligner::Node::Node(Node *parent, Kmer kmer, 
                      double prob, Type type) 
    :
      child_count_(0),
      kmer_(kmer),
      event_prob_(prob),
      type_(type),
      parent_(parent),
      children_(NULL) {

    #ifdef DEBUG_NODES
    full_length_ = parent->full_length_ + 1;
    #endif
    
    update_info();
}

//Creates invalid node
ForestAligner::Node::Node()
    : length_(0),
      stay_count_(0),
      skip_count_(0),
      ignore_count_(0),
      consec_stays_(0),
      child_count_(0),
      kmer_(0),
      event_prob_(0),
      seed_prob_(0),
      type_(Type::MATCH),
      parent_(NULL),
      children_(NULL) {}

//Copy constructor
ForestAligner::Node::Node(const Node &s)
    : length_(s.length_),
      stay_count_(s.stay_count_),
      skip_count_(s.skip_count_),
      ignore_count_(s.ignore_count_),
      consec_stays_(s.consec_stays_),
      child_count_(s.child_count_),
      kmer_(s.kmer_),
      event_prob_(s.event_prob_),
      seed_prob_(s.seed_prob_),
      type_(s.type_),
      parent_(s.parent_),
      children_(s.children_) {}
    
ForestAligner::Node::~Node() {
    if (children_ != NULL) {
        delete[] children_;
    }
}

bool ForestAligner::Node::is_valid() {
    return length_ > 0;
}

size_t ForestAligner::Node::seed_len() {
    return length_;
}

size_t ForestAligner::Node::match_len() {
    return seed_len() - stay_count_ - ignore_count_ + skip_count_;
}

double ForestAligner::Node::mean_prob() const {
    return seed_prob_ / length_;
}

double ForestAligner::Node::next_mean_prob(double next_prob) const {
    return (seed_prob_ + next_prob) / (length_ + 1);
}

//TODO: Really not a very good name
bool ForestAligner::Node::better_than(const Node *node) {
    return node->ignore_count_ > ignore_count_ ||
            (node->ignore_count_ == ignore_count_ && 
             node->mean_prob() < mean_prob());
}

bool ForestAligner::Node::should_report(const AlnParams &params) {
    bool r = type_ == Node::Type::MATCH && 
           stay_count_ <= params.max_stay_frac_ * length_ &&
           mean_prob() >= params.min_seed_pr_;
    //std::cout << stay_count_ << " stays " << length_ << " " << r << "\n";
    return r;
}

void ForestAligner::Node::invalidate(std::vector<Node *> *old_nodes, 
                                 bool delete_source = false) {

    Node *n = this, *p;

    while (n != NULL) {
        
        p = n->parent_;

        if (p != NULL || delete_source) {
            old_nodes->push_back(n);
        }

        n->child_count_ = 0;
        n->parent_ = NULL;
        n->length_ = n->stay_count_ = n->seed_prob_ = 0;

        if (p == NULL || p->remove_child(n)) {
            n = p;
        } else {
            break;
        }
    }
}

void ForestAligner::Node::update_info() {
    if (parent_ == NULL) {
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

    seed_prob_    = parent_->seed_prob_ + event_prob_;
    
    length_       = parent_->length_ + 1;
    stay_count_   = parent_->stay_count_;
    skip_count_   = parent_->skip_count_;
    ignore_count_ = parent_->ignore_count_;
    consec_stays_ = 0;

    switch(type_) {
        case Type::STAY:
        stay_count_++;
        consec_stays_ = parent_->consec_stays_ + 1;
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
    full_length_ = parent_->full_length_ + 1;
    #endif

    #ifdef DEBUG_PROB
    if (parent_->min_evt_prob_ < event_prob_) {
        min_evt_prob_ = parent_->min_evt_prob_;
    } else {
        min_evt_prob_ = event_prob_;
    }
    #endif
}


//Returns true if this node should be pruned
//If should_invalidate = false, will always return false
bool ForestAligner::Node::remove_child(Node *child) {

    if (child_count_ == 1) {
        child_count_ = 0;
        return true;
    }

    for (size_t c = 0; c < child_count_; c++) {
        if (children_[c] == child) {
            children_[c] = children_[--child_count_];
            break;
        }
    }

    return false;
}

size_t ForestAligner::Node::add_child(Node *child) {
    if (children_ == NULL) {
        children_ = new Node *[MAX_CHILDREN];
    }
    children_[child_count_++] = child;
    return child_count_;
}

void ForestAligner::Node::print() const {
    //std::cout << length_ << "\n";
              //<< full_length_ << "\n";
              //<< stay_count_ << "\t"
              ////<< skip_count_ << "\t"
              ////<< ignore_count_ << "\t"
              //<< (parents_.empty() ? Type::MATCH : parents_.front().second) << "\t"
              //<< parents_.size() << "\t"
              //<< children_.size() << "\t"
              //<< event_prob_ << "\t"
              //<< seed_prob_ << "\n";
}

ForestAligner::ForestAligner(const FMI &fmi, 
                     const AlnParams &ap,
                     const std::string &label)
    : fmi_(fmi),
      params_(ap),
      label_(label) {
    timer.reset();

    kmer_ranges_ = new Range[params_.model_.kmer_count()];
    for (Kmer k = 0; k < params_.model_.kmer_count(); k++) {
        Range r = fmi_.get_full_range(params_.model_.get_last_base(k));
        for (size_t i = params_.model_.kmer_len()-2; 
             i < params_.model_.kmer_len(); i--) {
            r = fmi_.get_neighbor(r, params_.model_.get_base(k, i));
        }
        kmer_ranges_[k] = r;
    }
}

ForestAligner::~ForestAligner() {
    delete[] kmer_ranges_;
    reset();
}

void ForestAligner::new_read(size_t read_len) {
    reset();
    cur_event_ = read_len;
}

void ForestAligner::reset() {
    for (auto p = prev_nodes_.begin(); p != prev_nodes_.end(); p++) {
        if (p->first->parent_ != NULL)
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

    timer.reset();
}

std::vector<Result> ForestAligner::add_event(double *kmer_probs, std::ostream &out) {

    //Update event index
    cur_event_--;

    Timer t;

    //Where this event's sources will be stored
    sources_.push_front(std::list<Node *>());
    
    double prob;

    t.reset();

    Range kmer_range;
    Kmer prev_kmer = 0;

    std::map<Range, Node *> prev_nodes_sort_;
    for (auto p = prev_nodes_.begin(); p != prev_nodes_.end(); p++) {
        prev_nodes_sort_[p->second] = p->first;
    }

    //Find neighbors of previous nodes
    //for (auto p = prev_nodes_.begin(); p != prev_nodes_.end(); p++) {
    for (auto p = prev_nodes_sort_.begin(); p != prev_nodes_sort_.end(); p++) {
        //Range prev_range = p->second;
        //Node *prev_node = p->first;
        Range prev_range = p->first;
        Node *prev_node = p->second;

        //std::cout << cur_event_ << "\t"
        //          << prev_range.start_ << "-" << prev_range.end_ << "\t"
        //          << prev_node->stay_count_ << "\t"
        //          << prev_node->length_ << "\t"
        //          << prev_node->mean_prob() << "\n";


        
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

        //Get probability for stay neighbor
        prev_kmer = prev_node->kmer_;
        prob = kmer_probs[prev_kmer];

        if (prev_node->consec_stays_ < params_.max_consec_stay_ && prob >= evpr_thresh) { //&& rank <= evpr_thresh ) { 
            neighbor_count += 1;
            //std::cout << prev_range.start_ << "-" << prev_range.end_ << "\t" << prob << "\n";

            Node next_node(prev_node, prev_kmer, prob, Node::Type::STAY);

            add_child(prev_range, next_node);
        }
        
        //Find next possible kmers
        auto neighbor_itr = params_.model_.get_neighbors(prev_kmer);
        std::list<Kmer> next_kmers;

        for (auto n = neighbor_itr.first; n != neighbor_itr.second; n++) {
            prob = kmer_probs[*n];
            if(prob >= evpr_thresh) { //&& prev_node->next_mean_prob(prob) >= params_.min_seed_pr_) {
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

            Base next_base = params_.model_.get_first_base(*next_kmer);
            Range next_range = fmi_.get_neighbor(prev_range, next_base);
            //std::cout << next_range.start_ << "-" << next_range.end_ << "\t" << prob << "\n";

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
    }

    //if (cur_event_ < 5000) {
    //    std::cout << timer.lap() << std::endl;
    //}

    //DEBUG(t.lap() << "\t");

    //Find sources
    for (Kmer kmer = 0; kmer < params_.model_.kmer_count(); kmer++) {
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

        if (prev_node->child_count_ == 0) {
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

    return r;
}

size_t ForestAligner::add_sources(const Range &range, const Node &node) {

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

ForestAligner::Node *ForestAligner::add_child(Range &range, Node &node) {
    
    if (!node.is_valid() || !range.is_valid()) {
        return NULL;
    }
    
    //Find closest node >= the node being added
    auto lb = next_nodes_.lower_bound(range);
    //std::cout << "a\n";

    //Node range hasn't been added yet
    if (lb == next_nodes_.end() || !lb->first.same_range(range)) {

        //Allocate memory for a new node
        Node *new_node;
        if (old_nodes_.empty()) {
            //std::cout << "b\n";
            new_node = new Node(node);
        } else {
            //std::cout << "c\n";
            *(old_nodes_.back()) = node;
            new_node = old_nodes_.back();
            old_nodes_.pop_back();
        }

        //Store it with it's range
        next_nodes_.insert(lb, std::pair<Range, Node *>(range, new_node));

        //Update the parent
        node.parent_->add_child(new_node);

        //Created one node
        return new_node;
    }
    
    //Node associated with same range
    Node *dup_node = lb->second;
    //std::cout << "d\n";



    if(node.better_than(dup_node)) {
        //std::cout << dup_node->length_ << "\t" << node.length_ << "\t" << dup_node->mean_prob() << "\t" << node.mean_prob() << "\n";
        //std::cout << "e\n";
        dup_node->parent_->remove_child(dup_node);
        node.parent_->add_child(dup_node);
        *dup_node = node;
        return dup_node;
    }

    //std::cout << "f\n";

    return dup_node;
}



std::vector<Result> ForestAligner::pop_seeds(std::ostream &out) { 

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

        for (size_t c = 0; c < aln_en->child_count_; c++) {

            #ifdef DEBUG_PATHS
            std::vector<Node::Type> p(params_.graph_elen_);
            p[0] = Node::Type::MATCH;
            p[1] = (*c)->type_;
            paths.push_back(p);
            #endif 

            aln_en->children_[c]->parent_ = NULL;
            next_sources.push_back(aln_en->children_[c]);
            to_visit.push_back(aln_en->children_[c]);

        }
        
        double prev_len = aln_en->seed_len();
        while (!to_visit.empty()) {
    
            Node *n = to_visit.front();
            to_visit.pop_front();

            #ifdef DEBUG_PATHS
            std::vector<Node::Type> p = paths.front();
            paths.pop_front();

            if (n->seed_len() > 2) {
                p[n->seed_len()-1] = n->type_;
            }
            #endif 

            if (n->seed_len() != prev_len) {
                prev_len = n->seed_len();
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

                        out << label_ << "\t";

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
                if (n->child_count_ > 0) { //Should always be true?
                    paths.push_back(p);

                    for (unsigned int i = 1; i < n->child_count_); i++) {
                        paths.push_back(std::vector<Node::Type>(p));
                    }   
                }
                #endif

                for (size_t c = 0; c < n->child_count_; c++) {
                    to_visit.push_back(n->children_[c]);
                }
            }

            n->length_--;
            
            n->update_info();
        }
        
        aln_ends.pop_front();
        aln_en->invalidate(&old_nodes_, true);
        //delete aln_en;
    }

    sources_.pop_back();

    return results;

}

void ForestAligner::print_graph(bool verbose) {
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
        } else if (n->parent_ == NULL) {
            source_counts[event]++;
        } else if (n->child_count_ > 1) {
            branch_counts[event]++;
        } else {
            linear_counts[event]++;
        }

        if (verbose) {
            std::cout << n << " " << event << " " 
                      << n->seed_len() << " ";
            
            if (n->parent_ == NULL) {
                std::cout << -1;
            } else {
                std::cout << n->parent_;
            }
        }

        for (size_t c = 0; c < n->child_count_; c++) {
            if (verbose) {
                std::cout << " " << n->children_[c];
            }
            to_visit.insert(std::pair<unsigned int, Node *>(event + 1, n->children_[c]));
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

    std::cout << "== total nodes: " 
              << (source_total + branch_total + linear_total + invalid_total) 
              << "\t==\n";

    std::cout << "== end ==\n";
}

