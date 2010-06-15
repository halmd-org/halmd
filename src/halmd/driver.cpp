/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <cmath>

#include <halmd/io/logger.hpp>
#include <halmd/driver.hpp>

using namespace boost;
using namespace std;

namespace halmd
{

/**
 * Assemble module options
 */
template <int dimension>
void driver<dimension>::options(po::options_description& desc)
{
    po::options_description group("Simulation");
    group.add_options()
        ("steps,s", po::value<uint64_t>()->default_value(10000),
         "number of simulation steps")
        ("time,t", po::value<double>(),
         "total simulation time")
        ;
    desc.add(group);
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void driver<dimension>::depends()
{
    modules::required<_Self, core_type>();
}

template <int dimension>
driver<dimension>::driver(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , core(modules::fetch<core_type>(vm))
{
    // parse options
    if (vm["steps"].defaulted() && !vm["time"].empty()) {
        time_ = vm["time"].as<double>();
        steps_ = static_cast<uint64_t>(round(time_ / core->integrator->timestep()));
    }
    else {
        steps_ = vm["steps"].as<uint64_t>();
        time_ = steps_ * core->integrator->timestep();
    }

    LOG("number of integration steps: " << steps_);
    LOG("integration time: " << time_);
}

/**
 * Run simulation
 */
template <int dimension>
void driver<dimension>::run()
{
    LOG("starting simulation");

    for (uint64_t i = 0; i < steps_; ++i) {
        core->mdstep();
    }

    LOG("finished simulation");
}

// explicit instantiation
template class driver<3>;
template class driver<2>;

template class module<driver<3> >;
template class module<driver<2> >;

} // namespace halmd
