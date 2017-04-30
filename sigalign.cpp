/* sigalign.cpp 
 * align oxford nanopore event-level signals directly to the reference without
 * basecalling */

#include <iostream>

#include "fast5.hpp"
#include "model_tools.hpp"

int main(int argc, char** argv) 
{
  if (argc < 2) { 
      std::cerr << "need a file. " << 
          "(usage: sigalign fast5_file_1 [fast5_file_2, ...]" << std::endl;
      return 1;
  }

  /* load the model into a vector:
   * model has size 4096*4. every group of four represents a 6mer, and is
   * ordered lexocographically from AAAAAA to TTTTTT. format:
   * [level_mean1, level_stdv1, sd_mean1, sd_stdv1, 
   *  level_mean2, level_stdv2, sd_mean2, sd_stdv2, 
   *  ...]
   */
  std::vector<double> model = load_model("r9_250bps.nucleotide.6mer.template.model");

  /* navigate through all the read files provided */
  for (int i = 1; i < argc; i++) {
      std::vector<Event> events;
      if (not fast5::File::is_valid_file(argv[i])) {
          std::cerr << "not a valid file. skipping ... " << std::endl;
          continue;
      }
      fast5::File f;
      try
      {
          f.open(argv[i]);
          assert(f.is_open());
          std::cout << "processing read " << argv[i] << std::endl;
          /* make sure that the fast5 contains an event sequence */
          if (f.have_eventdetection_events()) {
              events = f.get_eventdetection_events();
              /* to navigate through events:
              for (auto it = events.begin(); it < events.end(); it++) {
                  Event e = *it;
                  to get event mean: e.mean;
                  to get event stdv: e.stdv;
                  to get event start:  e.start;
              }
              */
          } else {
              std::cerr << "file " << argv[i] << " does not contain events. skipping..." << std::endl;
              continue;
          }
          ScaleParams scale = get_scale_params(model, events);
          /* TODO: start alignment, pass in scale, model (?), events */
      }
      catch (hdf5_tools::Exception& e)
      {
          std::cout << "hdf5 error: " << e.what() << std::endl;
      }
  }
}
