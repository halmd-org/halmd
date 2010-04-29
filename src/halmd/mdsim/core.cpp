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

#include <halmd/mdsim/core.hpp>
#include <halmd/util/logger.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

/**
 * Initialize simulation
 */
template <int dimension>
core<dimension>::core(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , force(module<mdsim::force<dimension> >::fetch(vm))
  , neighbor(module<mdsim::neighbor<dimension> >::fetch(vm))
  , sort(module<mdsim::sort<dimension> >::fetch(vm))
  , integrator(module<mdsim::integrator<dimension> >::fetch(vm))
  , position(module<mdsim::position<dimension> >::fetch(vm))
  , velocity(module<mdsim::velocity<dimension> >::fetch(vm))
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
    position->set();
    velocity->set();
    neighbor->update();
    force->compute();

    LOG("starting simulation");

    for (uint64_t i = 0; i < steps_; ++i) {
        integrator->integrate();
        if (neighbor->check()) {
            sort->order();
            neighbor->update();
        }
        force->compute();
        integrator->finalize();
    }

    LOG("finished simulation");
}

/**
 * Resolve module dependencies
 */
template <int dimension>
void core<dimension>::resolve(po::options const& vm)
{
    if (vm["dimension"].as<int>() != dimension) {
        throw inept_module();
    }
    module<mdsim::force<dimension> >::resolve(vm);
    module<mdsim::neighbor<dimension> >::resolve(vm);
    module<mdsim::sort<dimension> >::resolve(vm);
    module<mdsim::integrator<dimension> >::resolve(vm);
    module<mdsim::position<dimension> >::resolve(vm);
    module<mdsim::velocity<dimension> >::resolve(vm);
}

/**
 * Assemble module options
 */
template <int dimension>
po::options_description
core<dimension>::options()
{
    po::options_description desc;
    desc.add_options()
        ("steps,s", po::value<uint64_t>()->default_value(10000),
         "number of simulation steps")
        ("time,t", po::value<double>(),
         "total simulation time")
        ;
    return desc;
}

// explicit instantiation
template class core<3>;
template class core<2>;

} // namespace mdsim

template class module<mdsim::core<3> >;
template class module<mdsim::core<2> >;

} // namespace halmd
