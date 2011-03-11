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

#include <boost/foreach.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/sampler.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{

/**
 * Initialize simulation
 */
template <int dimension>
sampler<dimension>::sampler(
    shared_ptr<core_type> core
  , uint64_t steps
  , unsigned int statevars_interval
  , unsigned int trajectory_interval
)
  : core(core)
  , steps_(steps)
  , total_time_(steps_ * core->integrator->timestep())
  , statevars_interval_(statevars_interval)
  , trajectory_interval_(trajectory_interval)
{
    LOG("number of integration steps: " << steps_);
    LOG("integration time: " << total_time_);
}

/**
 * register module runtime accumulators
 */
template <int dimension>
void sampler<dimension>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.msv_output, "msv_output", "output of macroscopic state variables");
    profiler.register_runtime(runtime_.total, "total", "total simulation runtime");
}

/**
 * Run simulation
 */
template <int dimension>
void sampler<dimension>::run()
{
    {
        scoped_timer<timer> timer_(runtime_.total);

        LOG("setting up simulation box");
        prepare_observables(true);
        core->prepare();
        sample(true);

        LOG("starting simulation run");
        while (core->step_counter() < steps_) {
            // prepare observables in case of a sampling step
            prepare_observables(core->step_counter() + 1 == steps_);

            // perform complete MD integration step
            core->mdstep();

            // sample system state and properties,
            // force sampling after last integration step
            sample(core->step_counter() == steps_);
        }

        LOG("finished simulation run");
    }

    for_each(
        profiling_writers.begin()
      , profiling_writers.end()
      , bind(&profiling_writer_type::write, _1)
    );
}

/**
 * Sample system state and system properties
 */
template <int dimension>
void sampler<dimension>::sample(bool force)
{
    uint64_t step = core->step_counter();
    double time = step * core->integrator->timestep();
    bool is_sampling_step = false;

    if (!(step % statevars_interval_) || force) {
        BOOST_FOREACH (shared_ptr<observable_type> const& observable, observables) {
            observable->sample(time);
            is_sampling_step = true;
        }
        if (statevars_writer) {
            scoped_timer<timer> timer_(runtime_.msv_output);
            statevars_writer->write();
        }
    }

    // allow value 0 for trajectory_interval_
    if (((trajectory_interval_ && !(step % trajectory_interval_)) || force)
          && trajectory_writer) {
        trajectory_writer->append(time);
        is_sampling_step = true;
    }

    if (is_sampling_step)
        LOG_DEBUG("system state sampled at step " << step);
}

/**
 *  call prepare() for each observable if
 *  the current step is a sampling step
 */
template <int dimension>
void sampler<dimension>::prepare_observables(bool force)
{
    uint64_t step = core->step_counter();
    // increment step counter, function is invoked before core->mdstep()
    ++step;

    if (!(step % statevars_interval_) || force) {
        BOOST_FOREACH (shared_ptr<observable_type> const& observable, observables) {
            observable->prepare();
        }
    }
}

template <int dimension>
void sampler<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("sampler_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            class_<sampler, shared_ptr<runner>, runner>(class_name.c_str())
                .def(constructor<
                    shared_ptr<core_type>
                  , uint64_t
                  , unsigned int
                  , unsigned int
                >())
                .def("register_runtimes", &sampler::register_runtimes)
                .def_readwrite("observables", &sampler::observables)
                .def_readwrite("statevars_writer", &sampler::statevars_writer)
                .def_readwrite("phase_space", &sampler::phase_space)
                .def_readwrite("trajectory_writer", &sampler::trajectory_writer)
                .def_readwrite("profiling_writers", &sampler::profiling_writers)
                .property("steps", &sampler::steps)
                .property("total_time", &sampler::total_time)
                .property("trajectory_interval", &sampler::trajectory_interval)
                .property("statevars_interval", &sampler::statevars_interval)
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &sampler<3>::luaopen
    ]
    [
        &sampler<2>::luaopen
    ];
}

} // namespace

// explicit instantiation
template class sampler<3>;
template class sampler<2>;

} // namespace halmd
