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

#ifndef HALMD_OBSERVABLES_RUNTIME_ESTIMATE_HPP
#define HALMD_OBSERVABLES_RUNTIME_ESTIMATE_HPP

#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <stdint.h>
#include <string>
#include <sys/time.h>

#include <halmd/mdsim/clock.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {

/**
 * estimate remaining runtime until programme completion
 *
 * The wall clock time used for the completed simulation steps is
 * linearly extrapolated to the total number of simulation steps.
 */
class runtime_estimate
{
public:
    typedef halmd::mdsim::clock clock_type;
    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;
    typedef signal_type::connection connection_type;

    static void luaopen(lua_State* L);

    boost::shared_ptr<clock_type> clock;

    runtime_estimate(
        boost::shared_ptr<clock_type> clock
      , uint64_t total_steps
      , uint64_t step_start
    );

    // estimate remaining runtime and output to log file
    void sample() const;

    connection_type on_sample(slot_function_type const& slot)
    {
        return on_sample_.connect(slot);
    }

    //! returns estimate on remaining runtime in seconds based on number of completed simulation steps
    double value(uint64_t step) const;

    //! format time given in seconds, second argument specificies precision
    static std::string format_time(double time, unsigned int prec);

protected:
    /** simulation step counter at start */
    uint64_t step_start_;
    /** simulation step counter at finish */
    uint64_t step_stop_;
    /** wall clock time when simulation was started */
    timeval start_time_;

    signal_type on_sample_;
};

} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_RUNTIME_ESTIMATE_HPP */
