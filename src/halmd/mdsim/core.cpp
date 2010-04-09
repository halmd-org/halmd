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

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/host/forces/lj.hpp>
#include <halmd/mdsim/host/integrator/verlet.hpp>
#include <halmd/mdsim/host/neighbor.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/position/lattice.hpp>
#include <halmd/mdsim/host/random.hpp>
#include <halmd/mdsim/host/velocity/boltzmann.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
core<dimension, float_type>::core(options const& vm)
{
    shared_ptr<mdsim::particle<dimension, float_type> > particle;
    shared_ptr<mdsim::box<dimension, float_type> > box;
    shared_ptr<mdsim::force<dimension, float_type> > force;
    shared_ptr<mdsim::neighbor<dimension, float_type> > neighbor;
    shared_ptr<mdsim::random> random;
    shared_ptr<mdsim::position<dimension, float_type> > position;
    shared_ptr<mdsim::velocity<dimension, float_type> > velocity;
    shared_ptr<mdsim::integrator<dimension, float_type> > integrator;

    particle.reset(new mdsim::host::particle<dimension, float_type>(vm));
    box.reset(new mdsim::box<dimension, float_type>(particle, vm));
    force.reset(new mdsim::host::forces::lj<dimension, float_type>(particle, box, vm));
    neighbor.reset(new mdsim::host::neighbor<dimension, float_type>(particle, force, box, vm));
    random.reset(new mdsim::host::random(vm));
    position.reset(new mdsim::host::position::lattice<dimension, float_type>(particle, box, random, vm));
    velocity.reset(new mdsim::host::velocity::boltzmann<dimension, float_type>(particle, random, vm));
    integrator.reset(new mdsim::host::integrator::verlet<dimension, float_type>(particle, box, force, neighbor, vm));

    position->set();
    velocity->set();
    integrator->integrate(vm["steps"].as<uint64_t>());
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
