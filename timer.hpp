#ifndef TIMER_H
#define TIMER_H

//#include <time.h>
#include <chrono>

using namespace std::chrono;

class Timer {
    private:
        //clock_t start;
        high_resolution_clock::time_point start;

    public:
        inline Timer() {
            reset();
        }

        inline void reset() {	
            //start = clock();
            start = high_resolution_clock::now();
        }

        inline double get() {
            //return float( 1000.0 * (long double)(clock() - start) / 
            //              (long double)CLOCKS_PER_SEC );
            return (duration_cast< duration<double> > (high_resolution_clock::now() - start).count()) * 1000.0;
        }

        inline double lap() {
            double ret = get();
            reset();
            return ret;
        }
};

#endif
