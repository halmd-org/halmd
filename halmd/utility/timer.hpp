/*
 * Copyright © 2008-2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_UTILITY_TIMER_HPP
#define HALMD_UTILITY_TIMER_HPP

#include <boost/system/error_code.hpp>
#include <boost/system/system_error.hpp>
#include <ctime>
#include <limits>

namespace halmd {

/**
 * High-resolution timer
 *
 * This timer models boost::timer, but uses POSIX clock_gettime()
 * instead of C Standard Library clock() for better resolution.
 */
class timer
{
public:
    /**
     * start timer
     */
    timer();

    /**
     * restart timer
     */
    void restart();

    /**
     * returns elapsted time
     */
    double elapsed() const;

    static double elapsed_max();

    /**
     * returns resolution of timer
     */
    static double elapsed_min();

private:
    static void gettime(struct timespec& ts);

    static void getres(struct timespec& ts);

    struct timespec ts_;
};

inline timer::timer()
{
    timer::gettime(ts_);
}

inline void timer::restart()
{
    timer::gettime(ts_);
}

inline double timer::elapsed() const
{
    struct timespec ts;
    timer::gettime(ts);
    return (ts.tv_sec - ts_.tv_sec) + 1.e-9 * (ts.tv_nsec - ts_.tv_nsec);
}

inline double timer::elapsed_max()
{
    return std::numeric_limits<double>::max();
}

inline double timer::elapsed_min()
{
    struct timespec ts;
    timer::getres(ts);
    return ts.tv_sec + 1.e-9 * ts.tv_nsec;
}

inline void timer::gettime(struct timespec& ts)
{
    boost::system::error_code ec(
        clock_gettime(CLOCK_MONOTONIC, &ts)
      , boost::system::generic_category()
    );
    if (ec != boost::system::errc::success) {
        throw boost::system::system_error(ec);
    }
}

inline void timer::getres(struct timespec& ts)
{
    boost::system::error_code ec(
        clock_getres(CLOCK_MONOTONIC, &ts)
      , boost::system::generic_category()
    );
    if (ec != boost::system::errc::success) {
        throw boost::system::system_error(ec);
    }
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_TIMER_HPP */
