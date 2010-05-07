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

#include <algorithm>
#include <cmath>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/integrator/verlet.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace integrator
{

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::options(po::options_description& desc)
{
    desc.add_options()
        ("timestep,h", po::value<double>()->default_value(0.001),
         "simulation timestep")
        ;
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::resolve(po::options const& vm)
{
    module<particle_type>::required(vm);
    module<box_type>::required(vm);
}

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  , box(module<box_type>::fetch(vm))
{
    // parse options
    timestep_ = vm["timestep"].as<double>();

    LOG("integration timestep: " << timestep_);
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate()
{
    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type& v = particle->v[i] += particle->f[i] * static_cast<float_type>(timestep_half_);
        vector_type& r = particle->r[i] += v * static_cast<float_type>(timestep_);
        // enforce periodic boundary conditions
        // TODO: reduction is now to (-L/2, L/2) instead of (0, L) as before
        // check that this is OK
        particle->image[i] += box->reduce_periodic(r);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::finalize()
{
    for (size_t i = 0; i < particle->nbox; ++i) {
        particle->v[i] += particle->f[i] * static_cast<float_type>(timestep_half_);
    }
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class verlet<3, double>;
template class verlet<2, double>;
#else
template class verlet<3, float>;
template class verlet<2, float>;
#endif

}}} // namespace mdsim::host::integrator

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::integrator::verlet<3, double> >;
template class module<mdsim::host::integrator::verlet<2, double> >;
#else
template class module<mdsim::host::integrator::verlet<3, float> >;
template class module<mdsim::host::integrator::verlet<2, float> >;
#endif

} // namespace halmd
