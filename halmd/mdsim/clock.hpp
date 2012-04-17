/*
 * Copyright © 2011-2012 Peter Colberg
 * Copyright © 2011 Felix Höfling
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

#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <stdint.h> // uint64_t

#include <halmd/utility/signal.hpp>

namespace halmd {
namespace mdsim {

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
    /** difference between two simulation step counts */
    typedef int64_t step_difference_type;
    /** difference between two simulation times */
    typedef double time_difference_type;

    clock();

    /** advance clock by one step */
    void advance()
    {
        ++step_;
        // multiply instead of increment to avoid accumulated summation errors
        time_ = time_origin_ + (step_ - step_origin_) * timestep();
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

    /**
     * returns integration time step
     */
    time_type timestep() const;

    /**
     * set integration time step
     */
    void set_timestep(time_type timestep);

    /**
     * connect slot to set time step signal
     */
    connection on_set_timestep(boost::function<void (time_type)> const& slot)
    {
        return on_set_timestep_.connect(slot);
    }

private:
    /** step counter */
    step_type step_;
    /** simulation time */
    time_type time_;
    /** step counter at most recent call of set_timestep() */
    step_type step_origin_;
    /** simulation time at most recent call of set_timestep() */
    time_type time_origin_;
    /** integration time step */
    boost::optional<time_type> timestep_;
    /** set time step signal */
    signal<void (time_type)> on_set_timestep_;
};

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_CLOCK_HPP */
