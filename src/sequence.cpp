#include "sequence.hpp"

template <typename T>
std::string ValArray::to_string() const {
    std::stringstream ss;
    ss << "[";
    if (size() > 6) {
        for (size_t i = 0; i < 3; i++) {
            ss << at(i) << " ";
        }
        ss << "...";
        for (size_t i = size()-4; i < size(); i++) {
            ss << " " << at(i);
        }
    } else {
        for (size_t i = 0; i < size()-1; i++) {
            ss << at(i) << " ";
        }
        ss << at(size()-1);
    }
    ss << "]";
    return ss.str();
}
