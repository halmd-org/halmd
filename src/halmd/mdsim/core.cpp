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
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim
{

/**
 * Initialize simulation
 */
template <int dimension>
core<dimension>::core(options const& vm)
  : force(factory<mdsim::force<dimension> >::fetch(vm))
  , neighbor(factory<mdsim::neighbor<dimension> >::fetch(vm))
  , integrator(factory<mdsim::integrator<dimension> >::fetch(vm))
  , position(factory<mdsim::position<dimension> >::fetch(vm))
  , velocity(factory<mdsim::velocity<dimension> >::fetch(vm))
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
            neighbor->update();
        }
        force->compute();
        integrator->finalize();
    }

    LOG("finished simulation");
}

// explicit instantiation
template class core<3>;
template class core<2>;

}} // namespace halmd::mdsim
