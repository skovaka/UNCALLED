#ifndef TIMER_H
#define TIMER_H

#include <time.h>


class Timer {
    private:
        clock_t start;
    
    public:
        inline Timer() {
            start = clock();
        }

        inline void reset() {	
            start = clock();
        }

        inline float get() {
            return float( 1000.0 * (long double)(clock() - start) / 
                          (long double)CLOCKS_PER_SEC );
            //return (float) (clock() - start);
        }

        inline float lap() {
            float ret = get();
            reset();
            return ret;
        }
};

#endif
