#ifndef INCL_ARGPARSE
#define INCL_ARGPARSE

#include <string>
#include <map>

class ArgParse {
    public:
    enum Type {FLAG, INT, DOUBLE, STRING};

    private:
    typedef struct {
        Type type;
        std::string name, desc;
        void *value;
    } Arg;
    
    std::string desc_;
    std::map<char, Arg> args_;

    public:

    ArgParse(std::string desc);

    ~ArgParse();

    void parse_args(int argc, char **argv);

    bool add_flag(char c, std::string name, std::string desc);
    bool add_int(char c, std::string name, int def, std::string desc);
    bool add_double(char c, std::string name, double def, std::string desc);
    bool add_string(char c, std::string name, std::string def, std::string desc);

    std::string get_name(char c);
    std::string get_desc(char c);

    bool get_flag(char c);
    int get_int(char c);
    double get_double(char c);
    std::string get_string(char c);

    std::string get_param_str();


};


#endif
