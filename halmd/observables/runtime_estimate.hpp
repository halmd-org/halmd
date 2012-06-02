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

#include <boost/circular_buffer.hpp>
#include <lua.hpp>
#include <memory>
#include <string>

#include <halmd/mdsim/clock.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>

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
private:
    typedef halmd::signal<void ()> signal_type;

public:
    typedef halmd::mdsim::clock clock_type;
    typedef clock_type::step_type step_type;
    typedef signal_type::slot_function_type slot_function_type;

    static void luaopen(lua_State* L);

    runtime_estimate(
        std::shared_ptr<clock_type> clock
      , step_type total_steps
    );
    /** sample current real-time and simulation step */
    void sample();
    /** log estimate of remaining runtime */
    void estimate() const;
    /** format time given in seconds, second argument specificies precision */
    static std::string format_time(double time, unsigned int prec);

private:
    typedef std::pair<timer, step_type> timer_pair_type;

    /** simulation clock */
    std::shared_ptr<clock_type> clock_;
    /** simulation step counter at start */
    step_type step_start_;
    /** simulation step counter at finish */
    step_type step_stop_;
    /** pairs of timer and simulation step */
    boost::circular_buffer<timer_pair_type> timer_;
};

} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_RUNTIME_ESTIMATE_HPP */
