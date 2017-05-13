#include <fstream>
#include <sstream>
#include <numeric>
#include <tuple>
#include <string>
#include "model_tools.hpp"

std::pair<double, double> get_mean_stdv(std::vector<double> v) 
{
    /* see http://stackoverflow.com/a/7616783 */
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdv = std::sqrt(sq_sum / v.size());
    return std::make_pair(mean, stdv);
}

std::vector<double> load_model(std::string model_fname) 
{
    std::string line;
    double num;
    std::ifstream model_file(model_fname);
    std::vector<double> model;
    while(std::getline(model_file, line)) {
        std::stringstream parser(line);
        if (line[0] == '#'  || line.find("kmer") == 0) {
            continue;
        }
        std::string kmer;
        parser >> kmer;
        for (int i = 0; i < 4; i++) {
            parser >> num;
            model.push_back(num);
        }
    }
    return model;
}

/* uses Method of Moments of an event sequence to scale a given pore model */
std::vector<double> scale_model(std::vector<double> model, std::vector<Event> events)
{
    std::vector<double> scaled_model;
    std::vector<double> model_means;
    std::vector<double> event_means;
    std::vector<double> event_stdvs;
    for (int i = 0; i < 4096; i++) {
        model_means.push_back(model[4*i+0]);
    }
    for (auto it = events.begin(); it != events.end(); it++) {
        Event e = *it;
        event_means.push_back(e.mean);
        event_stdvs.push_back(e.stdv);
    }
    /* method of moments - get mean, stdv of event-level means */
    std::pair<double, double> events_mean_stdv = get_mean_stdv(event_means);
    std::pair<double, double> model_mean_stdv = get_mean_stdv(model_means);

    /* get scaling parameters */
    // double drift = 0.0;
    double var = 1.0;
    double scale_sd = 1.0;
    double var_sd = 1.0;
    /* TODO: calculate drift? */
    double shift = events_mean_stdv.first -  model_mean_stdv.first;
    double scale = events_mean_stdv.second / model_mean_stdv.second;

    /* scale model
     * see:
     * https://github.com/mateidavid/nanocall/blob/master/src/nanocall/Pore_Model.hpp#L126
     */
    for (int i = 0; i < 4096; i++) {
        double level_mean = model[4*i+0];
        double level_stdv = model[4*i+1];
        double stdv_mean = model[4*i+2];
        double stdv_stdv = model[4*i+3];
        double stdv_lambda = std::pow(stdv_mean, 3) / std::pow(stdv_stdv, 2);
        /* scale model level means */
        scaled_model.push_back(level_mean * scale + shift);
        /* scale model level stdvs */
        scaled_model.push_back(level_stdv * var);
        /* scale model stdv means */
        scaled_model.push_back(stdv_mean * scale_sd);
        /* scale model stdv lambdas */
        stdv_lambda = stdv_lambda * var_sd;
        /* scale model stdv stdvs */
        scaled_model.push_back(std::pow(std::pow(stdv_mean, 3) / stdv_lambda, 0.5));
    }
    return scaled_model;
}
    

ScaleParams get_scale_params(std::vector<double> model, std::vector<Event> events) 
{
    ScaleParams scale;
    std::vector<double> model_means;
    std::vector<double> event_means;
    std::vector<double> event_stdvs;
    for (int i = 0; i < 4096; i++) {
        model_means.push_back(model[4*i+0]);
    }
    for (auto it = events.begin(); it != events.end(); it++) {
        Event e = *it;
        event_means.push_back(e.mean);
        event_stdvs.push_back(e.stdv);
    }
    /* method of moments - get mean, stdv of event-level means */
    std::pair<double, double> events_mean_stdv = get_mean_stdv(event_means);
    std::pair<double, double> model_mean_stdv = get_mean_stdv(model_means);
    /* get scaling parameters */
    scale.shift = events_mean_stdv.first -  model_mean_stdv.first;
    scale.scale = events_mean_stdv.second / model_mean_stdv.second;
    scale.drift = 0.0;
    scale.var = 1.0;
    scale.scale_sd = 1.0;
    scale.var_sd = 1.0;
    /* TODO: calculate drift? */
    return scale;
}

std::vector<Event> simulate_read(std::vector<double> &model, 
                                 std::vector<mer_id> &ref, 
                                 int start, int end) {
    std::vector<Event> read(end - start);

    for (int i = start; i < end; i++) {
        read[i-start].mean = model[4*ref[i]];
        read[i-start].stdv = model[4*ref[i]+1];
        read[i-start].length = 4;
    }

    return read;
}
