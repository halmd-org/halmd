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

#include <halmd/mdsim/host/velocity/boltzmann.hpp>
#include <halmd/util/log.hpp>

namespace halmd { namespace mdsim { namespace host { namespace velocity
{

using namespace boost;
using namespace std;

template <int dimension, typename float_type>
boltzmann<dimension, float_type>::boltzmann(options const& vm)
    : _Base(vm)
    // dependency injection
    , particle(dynamic_pointer_cast<particle_type>(factory<mdsim::particle<dimension, float_type> >::fetch(vm)))
    , random(dynamic_pointer_cast<random_type>(factory<mdsim::random>::fetch(vm)))
{
    // parse options
    temp_ = vm["temperature"].as<float>();
}

/**
 * Initialise velocities from Maxwell-Boltzmann distribution
 */
template <int dimension, typename float_type>
void boltzmann<dimension, float_type>::set()
{
    float_type r, sigma = sqrt(temp_);
    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type& v = particle->v[i];
        if (dimension == 3) {
            if (i % 2) {
                v[0] = r;
                random->normal(v[1], v[2], sigma);
            }
            else {
                random->normal(v[0], v[1], sigma);
                random->normal(v[2], r, sigma);
            }
        }
        else {
            random->normal(v[0], v[1], sigma);
        }
    }
    LOG("assigned Maxwell-Boltzmann velocity distribution: T = " << temp_);
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class boltzmann<3, double>;
template class boltzmann<2, double>;
#else
template class boltzmann<3, float>;
template class boltzmann<2, float>;
#endif

}}}} // namespace halmd::mdsim::host::velocity
