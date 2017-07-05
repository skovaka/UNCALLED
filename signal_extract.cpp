/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>
#include <string>
#include <unordered_map>
#include <math.h>

#include "fast5.hpp"
#include "model_tools.hpp"
#include "timer.h"


int main(int argc, char** argv) 
{


    std::string model_fname = argv[1];

    /* load the model into a vector:
     * model has size 4096*4. every group of four represents a 6mer, and is
     * ordered lexocographically from AAAAAA to TTTTTT. format:
     * [level_mean1, level_stdv1, sd_mean1, sd_stdv1, 
     *  level_mean2, level_stdv2, sd_mean2, sd_stdv2, 
     *  ...]
     */
    std::vector<double> model = load_model(model_fname);

    /* navigate through all the read files provided */
    for (int fix = 2; fix < argc; fix++) {
        std::vector<Event> events;
        if (not fast5::File::is_valid_file(argv[fix])) {
            std::cerr << "<" << argv[fix] << "> is not a valid file. skipping... " << std::endl;
            continue;
        }
        fast5::File f;
        try
        {
            f.open(argv[fix]);
            assert(f.is_open());
            // std::cerr << "processing read " << argv[fix] << std::endl;
            //make sure that the fast5 contains an event sequence 
            if (f.have_eventdetection_events()) {
                events = f.get_eventdetection_events();
            } else {
                std::cerr << "file " << argv[fix] << " does not contain events. skipping..." << std::endl;
                continue;
            }
            ScaleParams scale = get_scale_params(model, events);

            for (int i = 0; i < events.size(); i++) {
                double model_mean = (events[i].mean / scale.scale) - scale.shift;
                double model_stdv = events[i].stdv / scale.var;
                //std::cout << model_mean << "\t" << model_stdv << std::endl;
                std::cout << i << " " << events[i].mean 
                          << "\t" << events[i].stdv 
                          << "\t" << (double) events[i].length << std::endl;
            }
 
        }
        catch (hdf5_tools::Exception& e)
        {
            std::cerr << "hdf5 error: " << e.what() << std::endl;
        }
    }
}
