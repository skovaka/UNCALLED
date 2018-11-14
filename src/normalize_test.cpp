#include <iostream>
#include "normalizer.hpp"

const std::string MODEL_DEF = "/home/skovaka/Dropbox/code/jhu/UNCALLED/src/uncalled/models/r94_5mers.txt";

bool open_fast5(const std::string &filename, fast5::File &file) {
    if (!fast5::File::is_valid_file(filename)) {
        std::cerr << "Error: '" << filename << "' is not a valid file \n";
    }

    try {
        file.open(filename);
        
        if (!file.is_open()) {  
            std::cerr << "Error: unable to open '" << filename << "'\n";
            return false;
        }

        return true;
        
    } catch (hdf5_tools::Exception& e) {
        std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
        return false;
    }

    return false;
}

int main(int argc, char **argv) {
    std::ifstream fast5_list(argv[1]);

    KmerModel model(MODEL_DEF, true);
    Normalizer win_norm(model, event_detection_defaults, 1000);
    EventDetector full_det(event_detection_defaults);

    std::string fast5_name;
    while (getline(fast5_list, fast5_name)) {

        //:std::cout << fast5_name << "\n";
        //

        fast5::File fast5;
        open_fast5(fast5_name, fast5);
        std::vector<float> samples = fast5.get_raw_samples();
        fast5.close();

        std::vector<Event> all_events = full_det.get_all_events(samples), half_events;

        for (u32 i = all_events.size()/2; i < all_events.size(); i++) {
            half_events.push_back(all_events[i]);
        }
       
        full_det.reset();
        Normalizer full_norm(model, event_detection_defaults, half_events.size());
        NormParams np = model.get_norm_params(half_events);
        std::cout << "Model norm:     " << np.scale << "\t" << np.shift << "\n";

        u32 i = 0;
        for (auto s : samples) {
            full_norm.add_sample(s);
            //if (win_norm.add_sample(s)) {
            //    np = win_norm.get_params();
            //    std::cout << i << "\t" << np.scale << "\t" << np.shift << "\n";
            //    i++;
            //}
        }
        np = full_norm.get_params();
        std::cout << "Full stream:    " << np.scale << "\t" << np.shift << "\n\n";
        win_norm.reset();
    }
}
