/* Realtime timers
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

#ifndef MDSIM_TIMER_HPP
#define MDSIM_TIMER_HPP

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <signal.h>
#include <sstream>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

namespace mdsim
{

/**
 * Realtime timer
 */
class real_timer
{
public:
    /**
     * start timer
     */
    void start()
    {
	gettimeofday(&m_start, NULL);
    }

    /**
     * stop timer
     */
    void stop()
    {
	gettimeofday(&m_stop, NULL);
    }

    /**
     * returns elapsed realtime in seconds
     */
    double elapsed()
    {
	timeval tv;
	timersub(&m_stop, &m_start, &tv);
	return tv.tv_sec + tv.tv_usec * 1.0E-6;
    }

    /**
     * format time given in seconds
     */
    static std::string format(double t)
    {
	std::ostringstream os;
	if (t < 60)
	    os << std::fixed << std::setprecision(1) << t << " s";
	else if (t < 3600)
	    os << std::fixed << std::setprecision(1) << (t / 60) << " min";
	else
	    os << std::fixed << std::setprecision(1) << (t / 3600) << " h";
	return os.str();
    }

    /**
     * output formatted time to stream
     */
    friend std::ostream& operator<<(std::ostream& os, real_timer& tm)
    {
	os << format(tm.elapsed());
	return os;
    }

private:
    timeval m_start, m_stop;
};


/**
 * Loop iterator with remaining runtime estimator
 */
template <typename T>
class iterator_timer
{
public:
    enum {
	INTERVAL = 15,
    };

public:
    iterator_timer(T const& value)
    {
	m_value = value;
	m_handler[0] = ::signal(SIGALRM, stop);
	// estimate may be triggered by user via USR1 signal
	m_handler[1] = ::signal(SIGUSR1, start);
	start();
    }

    ~iterator_timer()
    {
	alarm(0);
	::signal(SIGALRM, m_handler[0]);
	::signal(SIGUSR1, m_handler[1]);
    }

    static void start(int signum = 0)
    {
	m_timer.start();
	m_start = m_value;
	alarm(INTERVAL);
    }
    
    static void stop(int signum = 0)
    {
	m_timer.stop();
	double t = m_timer.elapsed() * (m_total - m_value) / (m_value - m_start);
	LOG("estimated remaining runtime: " << real_timer::format(t));
    }

    bool operator<(T const& value)
    {
	return (m_value < (m_total = value)) ? true : false;
    }

    T operator++(int)
    {
	return m_value++;
    }

    T operator++()
    {
	return ++m_value;
    }

    T& operator*()
    {
	return m_value;
    }

private:
    boost::array<sighandler_t, 2> m_handler;
    static T m_start, m_value, m_total;
    static real_timer m_timer;
};

template <typename T> T iterator_timer<T>::m_start;
template <typename T> T iterator_timer<T>::m_value;
template <typename T> T iterator_timer<T>::m_total;
template <typename T> real_timer iterator_timer<T>::m_timer;

} // namespace mdsim

#endif /* ! MDSIM_TIMER_HPP */
