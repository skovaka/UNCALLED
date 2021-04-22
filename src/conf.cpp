#include "conf.hpp"
        
#ifdef PYBIND
std::vector<std::string> Conf::_PARAM_GROUPS = {},
                         Conf::_GLOBAL_PARAMS = {};
#endif
