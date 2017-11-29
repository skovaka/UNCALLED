/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>

#include "event_detector.hpp"
//#include "kmer_model.hpp"
extern "C" {
    #include "scrappie_structures.h"
    #include "fast5_interface.h"
}

/*
typedef struct {
    uint64_t start;
    float length;
    float mean;
    float stdv;
    int pos;
    int state;
} event_t;
*/

int main(int argc, char** argv) {
    raw_table raw = read_raw(argv[1], true);

    EventDetector ed(event_detection_defaults);
    
    event_t e; 
    for (size_t i = 0; i < raw.n; i++) {
        e = ed.add_signal(raw.raw[i]);

        if (e.length > 0) {
            std::cout << e.length << "\t"
                      << e.mean << "\t" 
                      << e.stdv << "\t"
                      << e.start  << "\n"; 
        }
    }
}

