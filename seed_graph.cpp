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
#include "nano_fmi.hpp"
#include "boost/math/distributions/students_t.hpp"

//#define DEBUG(s)
#define DEBUG(s) do { std::cerr << s; } while (0)


//Source constructor
SeedGraph::Node::Node(mer_id kmer, double prob)
    : max_length_(1),
      stay_count_(0),
      kmer_(kmer),
      prob_(prob) {}

//Child constructor
SeedGraph::Node::Node(Node *parent, mer_id kmer, 
                            double prob, bool stay) 
    : max_length_(parent->max_length_ + 1),
      stay_count_(parent->stay_count_ + stay),
      kmer_(kmer),
      prob_(parent->prob_ + prob) {

    parents_.push_front(parent_ptr(parent, stay));
}

//Creates invalid node
SeedGraph::Node::Node()
    : max_length_(0),
      stay_count_(0),
      kmer_(0),
      prob_(0) {}

//Copy constructor
SeedGraph::Node::Node(const Node &s)
    : max_length_(s.max_length_),
      stay_count_(s.stay_count_),
      kmer_(s.kmer_),
      prob_(s.prob_),
      parents_(s.parents_),
      children_(s.children_) {}
    
SeedGraph::Node::~Node() {
    //invalidate();
}
bool SeedGraph::Node::is_valid() {
    return max_length_ > 0;
}

void SeedGraph::Node::invalidate(bool print) {

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
bool SeedGraph::Node::remove_child(Node *child, bool should_invalidate, bool print) {

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

SeedGraph::SeedGraph(const KmerModel &model,
                     const NanoFMI &fmi, 
                     const NormParams &norm_params, 
                     int seed_len, int read_len,
                     double event_prob,
                     double seed_prob,
                     double stay_prob,
                     double stay_frac)
    : model_(model),
      fmi_(fmi),
      norm_params_(norm_params),
      seed_length_(seed_len),
      cur_event_(read_len), 
      max_stays_( (int) (stay_frac * seed_len) ),
      event_prob_(event_prob),
      seed_prob_(seed_prob),
      stay_prob_(stay_prob) {
}


std::vector<Result> SeedGraph::add_event(Event e) {

    //bool stay = get_stay_prob(events[i], events[i+1]) >= STAY_THRESH;
    //bool stay = false;
    bool stay = true;
    double prob;

    
    event_kmer_probs_.push_front(new double[model_.kmer_count()]);
    double *kmer_probs = event_kmer_probs_.front();
    for (int kmer = 0; kmer < model_.kmer_count(); kmer++)
        kmer_probs[kmer] = model_.event_match_prob(e, kmer, norm_params_);
    
    sources_.push_front(std::list<Node *>());

    cur_event_--;

    //std::cout << "Creating children\n";
    //

    //DEBUG("ADDING EVENT " << cur_event_ << " (" << sources_.size() << ")\n");

    int source_count = 0;

    //auto prev_ranges = traversed_ranges.begin(); //Empty on first iteration
    for (auto p = prev_nodes.begin(); p != prev_nodes.end(); p++) {
        //Node *next_node = NULL;

        Range prev_range = p->first;
        Node *prev_node = p->second;


        int prev_kmer = prev_node->kmer_;
        prob = kmer_probs[prev_kmer];

        if (stay && prob >= event_prob_) {
            Node next_node(prev_node, prev_kmer, prob, true);
            Node *nn = add_child(prev_range, next_node);
        
            //if (nn->parents_.size() > 1) {
            //    //DEBUG("NEW NODE " << nn->max_length_ << " " << nn << " ");
            //    for (auto par = nn->parents_.begin(); par != nn->parents_.end(); par++)
            //        //DEBUG("\t" << par->first->max_length_);
            //    //DEBUG("\n");
            //}
        }

        //DEBUG(prev_node << "(" << prev_node->children_.size() << ")\n");
        
        auto neighbor_itr = model_.get_neighbors(prev_kmer);
        std::list<mer_id> next_kmers;

        for (auto n = neighbor_itr.first; n != neighbor_itr.second; n++) 
            if(kmer_probs[*n] >= event_prob_) 
                next_kmers.push_back(*n);


        std::list<Range> next_ranges = fmi_.get_neigbhors(prev_range, next_kmers);
            
        auto next_kmer = next_kmers.begin(); 
        auto next_range = next_ranges.begin();

        
        int child_count = 0;
        //std::cout << "Adding nodes\n";
        while (next_kmer != next_kmers.end()) {
            prob = kmer_probs[*next_kmer];

            Node next_node(prev_node, *next_kmer, prob, false);
            Node *nn = add_child(*next_range, next_node);
            child_count++;

            next_kmer++; //Maybe do this in-place up above
            next_range++;
        }

    }

    for (auto p = prev_nodes.begin(); p != prev_nodes.end(); p++) {
        Node *prev_node = p->second;

        if (prev_node->children_.empty()) {
            bool is_source = prev_node->parents_.empty();

            //DEBUG("x" << cur_event_ << "x " << prev_node->max_length_ << "\n");

            prev_node->invalidate(prev_node->max_length_ == 31);
            
            if (!is_source) {
                delete prev_node;
            } 
        }
    }


    for (mer_id kmer = 0; kmer < model_.kmer_count(); kmer++) {
        prob = kmer_probs[kmer];
        if (prob >= event_prob_) {
            Range next_range = fmi_.get_full_range(kmer);

            if (next_range.is_valid()) {
                Node next_node(kmer, prob);
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

int SeedGraph::add_sources(const Range &range, const Node &node) {

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

        next_nodes[range] = new Node(node);
        sources_.front().push_back(next_nodes[range]); //TODO: store node, no double query
        return 1;
    }

    std::list<Range> split_ranges;

    Range rr(range); //Right range

    //Node *node_copy;

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
        next_nodes[*r] = new Node(node);
        sources_.front().push_back(next_nodes[*r]);
    }

    return split_ranges.size();
}

SeedGraph::Node *SeedGraph::add_child(Range &range, Node &node) {
    
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
        Node *new_node = new Node(node);


        //DEBUG("ADDING " << new_node << " " << next_nodes.size() << "\n");

        //Store it with it's range
        auto r = next_nodes.insert(lb, std::pair<Range, Node *>(range, new_node));

        //DEBUG("NEW NODE " << new_node << " " << next_nodes.size() << "\n");

        //Update the parent
        new_parent.first->children_.push_front(new_node);

        //Created one node
        return new_node;
    }
    
    //Node associated with same range
    Node *dup_node = lb->second;
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
            dup_node->parents_.erase(dup_parent);

            //Copy seed info
            dup_node->max_length_ = node.max_length_;
            dup_node->stay_count_ = node.stay_count_;
            dup_node->prob_ = node.prob_;

            //Insert new parent in old parent's place
            dup_node->parents_.push_front(new_parent);
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
        } else if (new_parent.first->prob_ > dup_parent->first->prob_) {

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



std::vector<Result> SeedGraph::pop_seeds() { //Big result gathering loop

    //Is there a long enough seed?
    if (sources_.size() < seed_length_)
        return std::vector<Result>();

    //DEBUG("SEED POP\n");

    //TODO: Store prev_nodes like this?
    std::unordered_map<Node *, Range> node_ranges;
    for (auto n = prev_nodes.begin(); n != prev_nodes.end(); n++)
        node_ranges[n->second] = n->first;


    std::list<Node *> &aln_ends = sources_.back(),
                          &next_sources = *std::next(sources_.rbegin());

    while (!aln_ends.empty()) {
        Node *aln_en = aln_ends.front();

        if (!aln_en->is_valid()) {
            aln_ends.pop_front();
            delete aln_en;
            continue;
        }

        std::list<Node *> to_visit;
        for (auto c = aln_en->children_.begin(); c != aln_en->children_.end(); c++) {
            (*c)->parents_.clear();
            next_sources.push_back(*c);
            to_visit.push_back(*c);
        }
        
        double prev_len = aln_en->max_length_;
        auto kmer_probs = event_kmer_probs_.rbegin();

        while (!to_visit.empty()) {


            Node *n = to_visit.front();
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
                    Node *to_erase;

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
                if (n->stay_count_ <= max_stays_ && n->prob_ >= seed_prob_*seed_length_) {
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
                                  << r.prob << "\t" << n->stay_count_ << " =\n";
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

    event_kmer_probs_.pop_back();
    sources_.pop_back();

    return std::vector<Result>();

}

void SeedGraph::print_graph() {
    std::set< std::pair<int, Node *> > to_visit;
    auto event_sources = sources_.rbegin();

    std::cout << "== nodeptr \"event\" max_len #parents children... (" 
              << sources_.size() << " events) ==\n";

    for (auto s = event_sources->begin(); s != event_sources->end(); s++) 
        to_visit.insert( std::pair<int, Node *> (0, *s) );

    bool iter_sources = true;

    int source_count = 0;

    while (!to_visit.empty()) {

        if (iter_sources) {
            source_count++;
            event_sources++;
        }

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
        
        std::cout << n << " " << event << " " 
                  << n->max_length_ << " " 
                  << n->parents_.size();

        for (auto c = n->children_.begin(); c != n->children_.end(); c++) {
            std::cout << " " << *c;
            to_visit.insert(std::pair<int, Node *>(event + 1, *c));
        }

        std::cout << "\n";

    }

    std::cout << "== end ==\n";
}
