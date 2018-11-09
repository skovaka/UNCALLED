/* MIT License
 *
 * Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef TIMER_H
#define TIMER_H

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
