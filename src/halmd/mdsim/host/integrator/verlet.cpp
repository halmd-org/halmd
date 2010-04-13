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

#include <halmd/mdsim/host/integrator/verlet.hpp>
#include <halmd/util/log.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim { namespace host { namespace integrator
{

template <int dimension, typename float_type>
verlet<dimension, float_type>::verlet(options const& vm)
    : _Base(vm)
    // dependency injection
    , particle(dynamic_pointer_cast<particle_type>(module<mdsim::particle<dimension> >::fetch(vm)))
    , box(dynamic_pointer_cast<box_type>(module<mdsim::box<dimension> >::fetch(vm)))
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
        vector_type L = box->length();
        vector_type& image = particle->image[i];
        // enforce periodic boundary conditions
        for (size_t j = 0; j < dimension; ++j) {
            // assumes that particle position wraps at most once per time-step
            if (r[j] > L[j]) {
                r[j] -= L[j];
                image[j] += 1;
            }
            else if (r[j] < 0) {
                r[j] += L[j];
                image[j] -= 1;
            }
        }
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

}}}} // namespace halmd::mdsim::host::integrator
