#ifndef _INCL_MODEL_TOOLS
#define _INCL_MODEL_TOOLS

#include "fast5.hpp"

typedef struct ScaleParams {
    double shift = 0.0;
    double scale = 1.0;
    double drift = 0.0;
    double var = 1.0;
    double scale_sd = 1.0;
    double var_sd = 1.0;
} ScaleParams;

typedef fast5::EventDetection_Event Event;
std::vector<double> load_model(std::string model_fname);
std::vector<double> scale_model(std::vector<double> model, std::vector<Event> events);
ScaleParams get_scale_params(std::vector<double> model, std::vector<Event> events);
std::vector<Event> simulate_read(std::vector<double> &model, std::vector<mer_id> &ref, int start, int end);
#endif

