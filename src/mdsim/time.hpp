/* time.hpp
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MDSIM_TIME_HPP
#define MDSIM_TIME_HPP

#include <sys/time.h>
#include <time.h>
#include <signal.h>


namespace mdsim
{

/**
 * simple timer class for measuring elapsed real time
 */
class timer
{
private:
    /** total elapsed time */
    struct timeval tv_sum;
    /** last start time, in seconds and microseconds since the Epoch */
    struct timeval tv_start;

public:
    /**
     * timer constructor
     */
    timer()
    {
	timerclear(&tv_sum);
	timerclear(&tv_start);
    }

    /**
     * start timer
     */
    void start()
    {
	if (!timerisset(&tv_start)) {
	    gettimeofday(&tv_start, NULL);
	}
    }

    /**
     * stop timer
     */
    void stop()
    {
	struct timeval tv;

	if (timerisset(&tv_start)) {
	    gettimeofday(&tv, NULL);

	    timersub(&tv, &tv_start, &tv);
	    timeradd(&tv_sum, &tv, &tv_sum);

	    timerclear(&tv_start);
	}
    }

    /**
     * return total elapsed time in milliseconds
     */
    double elapsed()
    {
	return (tv_sum.tv_sec * 1.e3) + (tv_sum.tv_usec * 1.e-3);
    }

    /**
     * reset timer
     */
    void reset()
    {
	timerclear(&tv_sum);
	timerclear(&tv_start);
    }
};


/**
 * alarm clock
 */
struct alarm
{
    static void set(unsigned int seconds)
    {
	struct itimerval it;

	// clear default signal handler
	signal(SIGALRM, SIG_IGN);

	// current timer value
	it.it_value.tv_sec = seconds;
	it.it_value.tv_usec = 0;
	// next timer value
	it.it_interval.tv_sec = 0;
	it.it_interval.tv_usec = 0;

	setitimer(ITIMER_REAL, &it, NULL);
    }

    static bool get()
    {
	struct itimerval it;

	getitimer(ITIMER_REAL, &it);

	return !timerisset(&it.it_value);
    }
};

} // namespace mdsim

#endif /* ! MDSIM_TIME_HPP */
