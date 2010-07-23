/*
 * Copyright © 2008-2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#ifndef HALMD_UTILITY_TIMER_HPP
#define HALMD_UTILITY_TIMER_HPP

#include <boost/system/error_code.hpp>
#include <boost/system/system_error.hpp>
#include <ctime>

namespace halmd
{
namespace utility
{

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
    timer()
    {
        timer::gettime(ts_);
    }

    /**
     * restart timer
     */
    void restart()
    {
        timer::gettime(ts_);
    }

    /**
     * returns elapsted time
     */
    double elapsed() const
    {
        struct timespec ts;
        timer::gettime(ts);
        return (ts.tv_sec - ts_.tv_sec) + 1.e-9 * (ts.tv_nsec - ts_.tv_nsec);
    }

    /**
     * returns resolution of timer
     */
    static double resolution()
    {
        struct timespec ts;
        timer::getres(ts);
        return ts.tv_sec + 1.e-9 * ts.tv_nsec;
    }

private:
    static void gettime(struct timespec& ts)
    {
        boost::system::error_code ec(
            clock_gettime(CLOCK_MONOTONIC, &ts)
          , boost::system::get_posix_category()
        );
        if (ec != boost::system::errc::success) {
            throw boost::system::system_error(ec);
        }
    }

    static void getres(struct timespec& ts)
    {
        boost::system::error_code ec(
            clock_gettime(CLOCK_MONOTONIC, &ts)
          , boost::system::get_posix_category()
        );
        if (ec != boost::system::errc::success) {
            throw boost::system::system_error(ec);
        }
    }

    struct timespec ts_;
};

} // namespace utility

} // namespace halmd

#endif /* ! HALMD_UTILITY_TIMER_HPP */
