/*
 * Copyright © 2010-2011  Felix Höfling
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

#ifndef HALMD_SAMPLER_HPP
#define HALMD_SAMPLER_HPP

#include <lua.hpp>

#include <halmd/io/profiling/writer.hpp>
#include <halmd/io/statevars/writer.hpp>
#include <halmd/io/trajectory/writer.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/observables/observable.hpp>
#include <halmd/observables/trajectory.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd
{

template <int dimension>
class sampler
{
public:
    typedef mdsim::core<dimension> core_type;
    typedef observables::observable<dimension> observable_type;
    typedef io::statevars::writer<dimension> statevars_writer_type;
    typedef observables::trajectory<dimension> trajectory_type;
    typedef io::trajectory::writer<dimension> trajectory_writer_type;
    typedef io::profiling::writer profiling_writer_type;
    typedef utility::profiler profiler_type;

    static void luaopen(lua_State* L);

    sampler(
        boost::shared_ptr<core_type> core
      , uint64_t steps
      , unsigned int statevars_interval
      , unsigned int trajectory_interval
    );
    void run();
    void sample(bool force=false);
    void prepare_observables(bool force=false);
    void register_runtimes(profiler_type& profiler);

    boost::shared_ptr<core_type> core;
    std::vector<boost::shared_ptr<observable_type> > observables;
    boost::shared_ptr<statevars_writer_type> statevars_writer;
    boost::shared_ptr<trajectory_type> trajectory;
    boost::shared_ptr<trajectory_writer_type> trajectory_writer;
    std::vector<boost::shared_ptr<profiling_writer_type> > profiling_writers;

    // module runtime accumulator descriptions
    HALMD_PROFILING_TAG( msv_output_, "output of macroscopic state variables" );
    HALMD_PROFILING_TAG( total_, "total simulation runtime" );

    uint64_t steps()
    {
        return steps_;
    }

    double total_time()
    {
        return total_time_;
    }

    unsigned statevars_interval()
    {
        return statevars_interval_;
    }

    unsigned trajectory_interval()
    {
        return trajectory_interval_;
    }

private:
    /** total number of integration steps */
    uint64_t steps_;
    /** total integration time in MD units */
    double total_time_;
    // value from option --sampling-state-vars
    unsigned statevars_interval_;
    // value from option --sampling-trajectory
    unsigned trajectory_interval_;

    // list of profiling timers
    boost::fusion::map<
        boost::fusion::pair<msv_output_, accumulator<double> >
      , boost::fusion::pair<total_, accumulator<double> >
    > runtime_;
};

} // namespace halmd

#endif /* ! HALMD_SAMPLER_HPP */
