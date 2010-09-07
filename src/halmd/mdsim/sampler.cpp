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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/sampler.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace boost::fusion;
using namespace std;

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
        ("sampling-interval", po::value<unsigned>()->default_value(25),
         "sample system state every given number of integration steps")
        ;
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void sampler<dimension>::depends()
{
    modules::depends<_Self, thermodynamics_type>::optional();
    modules::depends<_Self, profiler_type>::required();
}

/**
 * Initialize simulation
 */
template <int dimension>
sampler<dimension>::sampler(modules::factory& factory, po::options const& vm)
  // dependency injection
  : thermodynamics(modules::fetch<thermodynamics_type>(factory, vm))
  , profiler(modules::fetch<profiler_type>(factory, vm))
  // store options
  , sampling_interval_(vm["sampling-interval"].as<unsigned>())
{
    // register module runtime accumulators
    profiler->register_map(runtime_);
}

/**
 * Sample system state and system properties
 */
template <int dimension>
void sampler<dimension>::sample(uint64_t step, double time)
{
    if (thermodynamics && !(step % sampling_interval_))
        thermodynamics->sample(time);
}

// explicit instantiation
template class sampler<3>;
template class sampler<2>;

} // namespace mdsim

template class module<mdsim::sampler<3> >;
template class module<mdsim::sampler<2> >;

} // namespace halmd
