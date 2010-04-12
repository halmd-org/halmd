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
#include <halmd/mdsim/host/forces/lj.hpp>
#include <halmd/mdsim/host/integrator/verlet.hpp>
#include <halmd/mdsim/host/neighbor.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/position/lattice.hpp>
#include <halmd/mdsim/host/random.hpp>
#include <halmd/mdsim/host/velocity/boltzmann.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim
{

/**
 * Initialize simulation
 */
template <int dimension, typename float_type>
core<dimension, float_type>::core(options const& vm)
  : particle(new mdsim::host::particle<dimension, float_type>(vm))
  , box(new mdsim::box<dimension, float_type>(particle, vm))
  , force(new mdsim::host::forces::lj<dimension, float_type>(particle, box, vm))
  , neighbor(new mdsim::host::neighbor<dimension, float_type>(particle, force, box, vm))
  , random(new mdsim::host::random(vm))
  , integrator(new mdsim::host::integrator::verlet<dimension, float_type>(particle, box, force, neighbor, vm))
  , position(new mdsim::host::position::lattice<dimension, float_type>(particle, box, random, vm))
  , velocity(new mdsim::host::velocity::boltzmann<dimension, float_type>(particle, random, vm))
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
template <int dimension, typename float_type>
void core<dimension, float_type>::run()
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
#ifndef USE_HOST_SINGLE_PRECISION
template class core<3, double>;
template class core<2, double>;
#else
template class core<3, float>;
template class core<2, float>;
#endif

}} // namespace halmd::mdsim
