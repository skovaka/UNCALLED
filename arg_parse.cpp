#include <iostream>
#include <unistd.h>
#include <iomanip>
#include <sstream>
#include "arg_parse.hpp"

ArgParse::ArgParse(std::string desc) {
    desc_ = desc;
}

ArgParse::~ArgParse() {

}

std::string ArgParse::get_param_str() {
    const std::string DELIM = ".";
    std::stringstream ss;

    for (auto a = args_.begin(); a != args_.end(); a++) {

        if (a->second.type != Type::STRING) {
            
            if (a != args_.begin()) {
                ss << DELIM;
            }

            ss << a->first;

            switch(a->second.type) {
                case Type::FLAG:
                ss << get_flag(a->first);
                break;

                case Type::INT:
                ss << get_int(a->first);
                break;
                
                case Type::DOUBLE:
                ss << std::fixed << std::setprecision(2) << get_double(a->first);
                break;

                default:
                std::cerr << "Error: bad type\n";
                break;
            }

            
        }
    }

    return ss.str();
}

void ArgParse::parse_args(int argc, char **argv) {

    std::string opt_str;

    for (auto a = args_.begin(); a != args_.end(); a++) {
        opt_str.push_back(a->first);

        if (a->second.type != Type::FLAG) {
            opt_str.push_back(':');
        }
    }

    char o;

    while ( (o = getopt(argc, argv, opt_str.c_str())) != -1 ) {

        if (args_.count(o) == 0) {
            std::cerr << "Error: unrecognized argument " << o << "\n";
            continue;
        }

        
        Arg &a = args_[o];
        switch(a.type) {
            case Type::FLAG:
            *( (bool *) a.value ) = true;
            break;

            case Type::INT:
            *( (int *) a.value ) = atoi(optarg);
            break;
            
            case Type::DOUBLE:
            *( (double *) a.value ) = atof(optarg);
            break;

            case Type::STRING:
            *( (std::string *) a.value ) = std::string(optarg);
            break;
            
            default:
            std::cerr << "Error: bad type\n";
            break;
        }
    }
}

bool ArgParse::add_flag(char c, std::string name, std::string desc="") {

    if (args_.count(c) > 0) {
        return false;
    }

    bool *val = new bool;
    *val = false;

    Arg a = {Type::FLAG, name, desc, (void *) val};
    args_[c] = a;

    return true;
}

bool ArgParse::add_int(char c, std::string name, 
                       int def, std::string desc="") {

    if (args_.count(c) > 0) {
        return false;
    }

    int *val = new int;
    *val = def;

    Arg a = {Type::INT, name, desc, (void *) val};
    args_[c] = a;

    return true;
}

bool ArgParse::add_double(char c, std::string name, 
                          double def, std::string desc="") {

    if (args_.count(c) > 0) {
        return false;
    }

    double *val = new double;
    *val = def;

    Arg a = {Type::DOUBLE, name, desc, (void *) val};
    args_[c] = a;

    return true;

}

bool ArgParse::add_string(char c, std::string name, 
                          std::string def, std::string desc="") {
    
    if (args_.count(c) > 0) {
        return false;
    }

    std::string *val = new std::string;
    *val = def;

    Arg a = {Type::STRING, name, desc, (void *) val};
    args_[c] = a;

    return true;
}

std::string ArgParse::get_name(char c) {
    return args_[c].name;
}

std::string ArgParse::get_desc(char c) {
    return args_[c].desc;
}

bool ArgParse::get_flag(char c) {
    Arg &a = args_[c];
    return *( (bool *) a.value );
}

int ArgParse::get_int(char c) {
    Arg &a = args_[c];
    return *( (int *) a.value );
}

double ArgParse::get_double(char c) {
    Arg &a = args_[c];
    return *( (double *) a.value );
}

std::string ArgParse::get_string(char c) {
    Arg &a = args_[c];
    return *( (std::string *) a.value );
}

//int main(int argc, char** argv) {
//    ArgParse ap("This is a test");
//    ap.add_flag('f', "flag", "This is a flag");
//    ap.add_int('i', "int", 0, "This is an integer argument");
//    ap.add_double('d', "double", 0, "This is a double argument");
//    ap.add_string('s', "string", "", "This is a string argument");
//
//    ap.parse_args(argc, argv);
//
//    std::cout << ap.get_name('f') << ": " << ap.get_flag('f') << "\n";
//    std::cout << ap.get_name('i') << ": " << ap.get_int('i') << "\n";
//    std::cout << ap.get_name('d') << ": " << ap.get_double('d') << "\n";
//    std::cout << ap.get_name('s') << ": " << ap.get_string('s') << "\n";
//
//}
