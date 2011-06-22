/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#ifndef HALMD_MDSIM_CLOCK_HPP
#define HALMD_MDSIM_CLOCK_HPP

#include <stdint.h> // uint64_t

namespace halmd
{
namespace mdsim
{

/**
 * The clock module contains the current simulation step and time.
 */
class clock
{
public:
    /** simulation step counter type */
    typedef uint64_t step_type;
    /** simulation time type */
    typedef double time_type;

    clock(time_type timestep);

    /** advance clock by one step */
    void advance()
    {
        ++step_;
        // multiply instead of increment to avoid accumulated summation errors
        time_ = step_ * timestep_;
    }

    /** MD step counter */
    step_type step() const
    {
        return step_;
    }

    /** MD time */
    time_type time() const
    {
        return time_;
    }

    /** MD timestep */
    time_type timestep() const
    {
        return timestep_;
    }

private:
    /** step counter */
    step_type step_;
    /** simulation time */
    time_type time_;
    /** timestep */
    time_type timestep_;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_CLOCK_HPP */
