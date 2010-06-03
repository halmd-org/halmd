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
#include <halmd/mdsim/core.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Assemble module options
 */
template <int dimension>
void core<dimension>::options(po::options_description& desc)
{
    desc.add_options()
        ("dimension", po::value<int>()->default_value(3),
         "positional coordinates dimension")
        ("steps,s", po::value<uint64_t>()->default_value(10000),
         "number of simulation steps")
        ("time,t", po::value<double>(),
         "total simulation time")
        ;
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void core<dimension>::resolve(po::options const& vm)
{
    if (vm["dimension"].as<int>() != dimension) {
        throw unsuitable_module("mismatching option 'dimension'");
    }
    module<force_type>::required(vm);
    module<neighbor_type>::required(vm);
    module<sort_type>::optional(vm);
    module<integrator_type>::required(vm);
    module<position_type>::required(vm);
    module<velocity_type>::required(vm);
}

/**
 * Initialize simulation
 */
template <int dimension>
core<dimension>::core(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , force(module<force_type>::fetch(vm))
  , neighbor(module<neighbor_type>::fetch(vm))
  , sort(module<sort_type>::fetch(vm))
  , integrator(module<integrator_type>::fetch(vm))
  , position(module<position_type>::fetch(vm))
  , velocity(module<velocity_type>::fetch(vm))
{
    // parse options
    if (vm["steps"].defaulted() && !vm["time"].empty()) {
        time_ = vm["time"].as<double>();
        steps_ = static_cast<uint64_t>(round(time_ / integrator->timestep()));
    }
    else {
        steps_ = vm["steps"].as<uint64_t>();
        time_ = steps_ * integrator->timestep();
    }

    LOG("number of integration steps: " << steps_);
    LOG("integration time: " << time_);
}

/**
 * Run simulation
 */
template <int dimension>
void core<dimension>::run()
{
    init();

    LOG("starting simulation");

    for (uint64_t i = 0; i < steps_; ++i) {
        mdstep();
    }

    LOG("finished simulation");
}

/**
 * Initialise simulation
 */
template <int dimension>
inline void core<dimension>::init()
{
    position->set();
    velocity->set();
    neighbor->update();
    force->compute();
}

/**
 * Perform a single MD integration step
 */
template <int dimension>
inline void core<dimension>::mdstep()
{
    integrator->integrate();
    if (neighbor->check()) {
        if (sort) {
            sort->order();
        }
        neighbor->update();
    }
    force->compute();
    integrator->finalize();
}

// explicit instantiation
template class core<3>;
template class core<2>;

} // namespace mdsim

template class module<mdsim::core<3> >;
template class module<mdsim::core<2> >;

} // namespace halmd
