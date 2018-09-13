#include <iostream>
#include <string>
#include "fast5.hpp"
#include "basepairs.hpp"

void dump_bcevents(fast5::File &file) {
    std::vector<fast5::Basecall_Event> bcevents = file.get_basecall_events(0);
    for (auto e : bcevents) {
        std::cout << e.mean << "\t";
                  //<< e.stdv << "\t";

        //if (e.length < 1) {
        //    std::cout << (int) (4000*e.length) << "\t";
        //} else {
        //    std::cout << e.length << "\t";
        //}
        
        for (auto c : e.model_state) {
            if (BASE_BYTES[c] < 4) std::cout << c;
        }
        //std::cout << "\t" << e.move << std::endl;
        std::cout << std::endl;
        //std::cout.flush();
    }

}

void dump_raw_info(fast5::File &file) {
    auto params = file.get_raw_samples_params();
    std::cout << params.read_number << "\t"
              << params.start_mux << "\t"
              << params.start_time << "\t"
              << params.duration << std::endl;
}

int main(int argc, char** argv) {
    std::ifstream reads_in(argv[1]);
    std::string fast5_name;

    while (getline(reads_in, fast5_name)) {
        if (!fast5::File::is_valid_file(fast5_name)) {
            std::cerr << "Error: '" << fast5_name << "' is not a valid file \n";
        }

        try {
            fast5::File file;
            file.open(fast5_name);
            
            if (!file.is_open()) {  
                std::cerr << "Error: unable to open '" 
                          << fast5_name << "'\n";
                continue;
            }

            std::cout << "== " << fast5_name << " ==" << std::endl;
            dump_raw_info(file);
            dump_bcevents(file);
            std::cout << "== END ==" << std::endl;

        } catch (hdf5_tools::Exception& e) {
            std::cerr << "Error: hdf5 exception '" << e.what() << "'\n";
        }
    }
}
