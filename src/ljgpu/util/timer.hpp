/* Realtime timers
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_UTIL_TIMER_HPP
#define LJGPU_UTIL_TIMER_HPP

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <signal.h>
#include <sstream>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

namespace ljgpu
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
	else if (t < 86400)
	    os << std::fixed << std::setprecision(1) << (t / 3600) << " h";
	else
	    os << std::fixed << std::setprecision(1) << (t / 86400) << " d";
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
 * Iterator with remaining realtime estimator
 */
template <typename T, time_t I = 15, time_t W = 60>
class count_timer
{
public:
    /**
     * initialize iterator variable and schedule timer
     */
    count_timer(T const& i = T()) :
	m_start(time(NULL) + W),
	m_stop(-1),
	m_count(i),
	m_time(0) {}

    /**
     * schedule timer after given time interval in seconds
     */
    void set(time_t wait = W)
    {
	m_start = time(NULL) + wait;
	m_stop = -1;
	m_time = 0;
    }

    /**
     * reset estimated remaining time
     */
    void clear()
    {
	m_start = -1;
	m_stop = -1;
	m_time = 0;
    }

    /**
     * returns estimated remaining time
     */
    double const& estimate()
    {
	return m_time;
    }

    operator T& ()
    {
	return m_count;
    }

    operator T const& () const
    {
	return m_count;
    }

    template <typename U>
    bool operator<(U const& value)
    {
	handle_timer(value);
	return (m_count < value);
    }

    template <typename U>
    bool operator>(U const& value)
    {
	handle_timer(value);
	return (m_count > value);
    }

    template <typename U>
    bool operator<=(U const& value)
    {
	handle_timer(value + 1);
	return (m_count <= value);
    }

    template <typename U>
    bool operator>=(U const& value)
    {
	handle_timer(value - 1);
	return (m_count >= value);
    }

    template <typename U>
    friend bool operator<(U const& value, count_timer<T>& timer)
    {
	return (timer > value);
    }

    template <typename U>
    friend bool operator>(U const& value, count_timer<T>& timer)
    {
	return (timer < value);
    }

    template <typename U>
    friend bool operator<=(U const& value, count_timer<T>& timer)
    {
	return (timer >= value);
    }

    template <typename U>
    friend bool operator>=(U const& value, count_timer<T>& timer)
    {
	return (timer <= value);
    }

private:
    void handle_timer(T const& limit)
    {
	// time() is about an order of magnitude faster than gettimeofday()
	const time_t t = time(NULL);

	if (m_start != -1 && t >= m_start) {
	    m_timer.start();
	    m_stop = t + I;
	    m_start = -1;
	    m_base = m_count;
	}
	else if (m_stop != -1 && t >= m_stop) {
	    m_timer.stop();
	    m_stop = -1;
	    m_time = m_timer.elapsed() * (limit - m_count) / (m_count - m_base);
	}
    }

private:
    time_t m_start, m_stop;
    real_timer m_timer;
    T m_base, m_count;
    double m_time;
};

/**
 * returns time interval in seconds between two timespecs
 */
inline double operator-(timespec const& ts1, timespec const& ts2)
{
    return (ts1.tv_sec - ts2.tv_sec) + (ts1.tv_nsec - ts2.tv_nsec) * 1.E-9;
}

} // namespace ljgpu

#endif /* ! LJGPU_UTIL_TIMER_HPP */
