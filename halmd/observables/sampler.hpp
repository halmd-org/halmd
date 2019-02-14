/*
 * Copyright © 2010-2012 Peter Colberg
 * Copyright © 2010-2011 Felix Höfling
 * Copyright © 2013      Nicolas Höft
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

#ifndef HALMD_OBSERVABLES_SAMPLER_HPP
#define HALMD_OBSERVABLES_SAMPLER_HPP

#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {
namespace observables {

/**
 * Sampler to run Molecular Dynamics simulation
 */
class sampler
{
public:
    typedef mdsim::clock clock_type;
    typedef clock_type::step_type step_type;
    typedef mdsim::core core_type;

    sampler(
        std::shared_ptr<clock_type> clock
      , std::shared_ptr<core_type> core
    );

    /**
     * Sample current state
     */
    void sample();

    /**
     * Run simulation for given number of steps
     */
    void run(step_type steps);

    /**
     * Connect slot to signal emitted before MD integration step
     */
    connection on_prepare(std::function<void ()> const& slot, step_type interval, step_type start);

    /**
     * Connect slot to signal emitted to sample current state
     */
    connection on_sample(std::function<void ()> const& slot, step_type interval, step_type start);

    /**
     * Connect slot to signal emitted before starting simulation run
     */
    connection on_start(std::function<void ()> const& slot);

    /**
     * Connect slot to signal emitted after finishing simulation run
     */
    connection on_finish(std::function<void ()> const& slot);

    /**
     * Bind class to Lua
     */
    static void luaopen(lua_State* L);

private:
    /** simulation clock */
    std::shared_ptr<clock_type> clock_;
    /** simulation core */
    std::shared_ptr<core_type> core_;
    /** signal emitted before MD integration step */
    signal<void ()> on_prepare_;
    /** signal emitted after MD integration step */
    signal<void ()> on_sample_;
    /** signal emitted before starting simulation run */
    signal<void ()> on_start_;
    /** signal emitted after finishing simulation run */
    signal<void ()> on_finish_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type total;
        accumulator_type prepare;
        accumulator_type sample;
        accumulator_type start;
        accumulator_type finish;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLER_HPP */
