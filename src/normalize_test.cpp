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
    fast5::File fast5;
    open_fast5(argv[1], fast5);
    std::vector<float> samples = fast5.get_raw_samples();

    KmerModel model(MODEL_DEF, true);
    Normalizer norm(model, event_detection_defaults, 4000);

    for (auto s : samples) {
        if (norm.add_sample(s))
            std::cout << norm.pop_event() << " n\n";
    }

    std::cout << "It works " << samples.size() << "\n";
}
