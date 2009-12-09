/* Exception class
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_UTIL_EXCEPTION_HPP
#define HALMD_UTIL_EXCEPTION_HPP

#include <exception>

namespace halmd {

class exception : public std::exception
{
public:
    exception(char const* str): str_(str) {}

    char const* what() const throw()
    {
        return str_;
    }

private:
    char const* str_;
};

/**
 * exit codes
 */
enum {
    // no errors
    HALMD_EXIT_SUCCESS = 0,
    // aborted simulation on SIGTERM or SIGINT
    HALMD_EXIT_TERM = 1,
    // caught standard exception
    HALMD_EXIT_EXCEPTION = 2,
    // caught CUDA error
    HALMD_EXIT_CUDA_ERROR = 3,
    // potential energy diverged
    HALMD_EXIT_POTENTIAL_ENERGY_DIVERGENCE = 4,
};

} // namespace halmd

#endif /* ! HALMD_UTIL_EXCEPTION_HPP */
