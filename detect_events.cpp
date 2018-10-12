/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <iomanip>
#include <string>

#include "event_detector.hpp"

int main(int argc, char** argv) {

    std::string filename = std::string(argv[1]);

    detector_param params = event_detection_defaults;
    params.threshold1 = 1.4;
    params.threshold2 = 1.1;
    EventDetector ed(params, 30, 150);

    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        fast5::File file;
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" 
                      << filename << "'\n";
        }

        std::vector<Event> events = ed.get_all_events(file.get_raw_samples());

        std::cout << std::setprecision(4);

        for (size_t i = 0; i < events.size(); i++) {
            Event e = events[i];
            std::cout << e.mean << "\t" 
                      << e.stdv << "\t"
                      << e.length << "\t"
                      << e.start  << "\n"; 
        }
        

    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
    }
}

