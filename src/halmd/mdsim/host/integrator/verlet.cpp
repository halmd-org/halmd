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
verlet<dimension, float_type>::verlet(particle_ptr particle, box_ptr box, force_ptr force, neighbor_ptr neighbor, options const& vm)
    // dependency injection
    : particle(static_pointer_cast<particle_type>(particle))
    , box(static_pointer_cast<box_type>(box))
    , force(static_pointer_cast<force_type>(force))
    , neighbor(static_pointer_cast<neighbor_type>(neighbor))
{
    // parse options
    timestep_ = vm["timestep"].as<double>();

    LOG("integration timestep: " << timestep_);
}

template <int dimension, typename float_type>
void verlet<dimension, float_type>::integrate(uint64_t steps)
{
    neighbor->update();
    force->compute();
    for (uint64_t i = 0; i < steps; ++i) {
        pre_force();
        if (neighbor->check()) {
            neighbor->update();
        }
        force->compute();
        post_force();
    }
}

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type>
void verlet<dimension, float_type>::pre_force()
{
    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type& v = particle->v[i] += particle->f[i] * timestep_half_;
        vector_type& r = particle->r[i] += v * timestep_;
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
void verlet<dimension, float_type>::post_force()
{
    for (size_t i = 0; i < particle->nbox; ++i) {
        particle->v[i] += particle->f[i] * timestep_half_;
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
