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
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/integrator/verlet.hpp>
#include <halmd/utility/lua.hpp>
#include <halmd/utility/module.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace integrator
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::depends()
{
    modules::depends<_Self, particle_type>::required();
    modules::depends<_Self, box_type>::required();
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::select(po::variables_map const& vm)
{
    if (vm["integrator"].as<std::string>() != "verlet") {
        throw unsuitable_module("mismatching option integrator");
    }
}

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(modules::factory& factory, po::variables_map const& vm)
  : _Base(factory, vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(factory, vm))
  , box(modules::fetch<box_type>(factory, vm))
  // set parameters
  , timestep_half_(0.5 * timestep_)
{
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
