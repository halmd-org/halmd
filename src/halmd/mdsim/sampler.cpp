/*
 * Copyright © 2010  Felix Höfling
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
#include <halmd/mdsim/sampler.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost::fusion;

namespace halmd
{
namespace mdsim
{

/**
 * Assemble module options
 */
template <int dimension>
void sampler<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("sampling-state-vars", po::value<unsigned>()->default_value(25),
         "sample macroscopic state variables every given number of integration steps")
        ;
    desc.add_options()
        ("sampling-trajectory", po::value<unsigned>()->default_value(0),
         "sample trajectory every given number of integration steps")
        ;
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void sampler<dimension>::depends()
{
    modules::depends<_Self, core_type>::required();
    modules::depends<_Self, observable_type>::optional();
    modules::depends<_Self, statevars_writer_type>::optional();
    modules::depends<_Self, trajectory_writer_type>::required();
}

/**
 * Initialize simulation
 */
template <int dimension>
sampler<dimension>::sampler(modules::factory& factory, po::options const& vm)
  // dependency injection
  : core(modules::fetch<core_type>(factory, vm))
  , observables(modules::fetch<observable_type>(factory, vm))
  , statevars_writer(modules::fetch<statevars_writer_type>(factory, vm))
  , trajectory_writer(modules::fetch<trajectory_writer_type>(factory, vm))
  // store options
  , statevars_interval_(vm["sampling-state-vars"].as<unsigned>())
  , trajectory_interval_(vm["sampling-trajectory"].as<unsigned>())
{
    /*@{ FIXME remove pre-Lua hack */
    shared_ptr<profiler_type> profiler(modules::fetch<profiler_type>(factory, vm));
    register_runtimes(*profiler);
    /*@}*/

    BOOST_FOREACH (shared_ptr<observable_type> const& observable, observables) {
        observable->register_statevars(*statevars_writer);
    }
}

/**
 * register module runtime accumulators
 */
template <int dimension>
void sampler<dimension>::register_runtimes(profiler_type& profiler)
{
    profiler.register_map(runtime_);
}

/**
 * Sample system state and system properties
 */
template <int dimension>
void sampler<dimension>::sample(bool force)
{
    uint64_t step = core->step_counter();
    bool is_sampling_step = false;

    if (!(step % statevars_interval_) || force) {
        BOOST_FOREACH (shared_ptr<observable_type> const& observable, observables) {
            observable->sample(core->time());
            is_sampling_step = true;
        }
        if (statevars_writer) {
            scoped_timer<timer> timer_(at_key<msv_output_>(runtime_));
            statevars_writer->write();
        }
    }

    // allow value 0 for trajectory_interval_
    if ((trajectory_interval_ && !(step % trajectory_interval_) || force)
          && trajectory_writer) {
        trajectory_writer->append();
        is_sampling_step = true;
    }

    if (is_sampling_step)
        LOG_DEBUG("system state sampled at step " << step);
}

// explicit instantiation
template class sampler<3>;
template class sampler<2>;

} // namespace mdsim

template class module<mdsim::sampler<3> >;
template class module<mdsim::sampler<2> >;

} // namespace halmd
