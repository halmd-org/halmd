/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#ifndef HALMD_OBSERVABLES_SAMPLER_HPP
#define HALMD_OBSERVABLES_SAMPLER_HPP

#include <stdint.h>
#include <lua.hpp>
#include <utility> // pair

#include <halmd/mdsim/core.hpp>
#include <halmd/runner.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd
{
namespace observables
{

/**
 * Sampler to run Molecular Dynamics simulation
 */
template <int dimension>
class sampler
  : public runner
{
public:
    typedef mdsim::core<dimension> core_type;
    typedef utility::profiler profiler_type;
    typedef halmd::signal<void (uint64_t)> signal_type;
    typedef typename signal_type::slot_function_type slot_function_type;

    struct runtime
    {
        typedef typename profiler_type::accumulator_type accumulator_type;
        accumulator_type total;
    };

    sampler(boost::shared_ptr<core_type> core, uint64_t steps);
    virtual void run();
    void register_runtimes(profiler_type& profiler);
    void on_start(slot_function_type const& slot);
    void on_prepare(slot_function_type const& slot, uint64_t interval);
    void on_sample(slot_function_type const& slot, uint64_t interval);
    void on_finish(slot_function_type const& slot);
    static void luaopen(lua_State* L);

    uint64_t steps()
    {
        return steps_;
    }

    double total_time()
    {
        return total_time_;
    }

private:
    void prepare(slot_function_type const& slot, uint64_t interval, uint64_t step) const;
    void sample(slot_function_type const& slot, uint64_t interval, uint64_t step) const;

    /** Molecular Dynamics simulation core */
    boost::shared_ptr<core_type> core_;
    /** total number of integration steps */
    uint64_t steps_;
    /** total integration time in MD units */
    double total_time_;
    /** profiling runtime accumulators */
    runtime runtime_;
    /** signal emitted before starting simulation run */
    signal_type on_start_;
    /** signal emitted before MD integration step */
    signal_type on_prepare_;
    /** signal emitted after MD integration step */
    signal_type on_sample_;
    /** signal emitted after finishing simulation run */
    signal_type on_finish_;
};

} // namespace observables

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_SAMPLER_HPP */
