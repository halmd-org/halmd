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

#include <halmd/io/logger.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/particle.hpp>
#endif
#include <halmd/mdsim/samples/host/trajectory.hpp>

namespace halmd
{
namespace mdsim { namespace samples { namespace host
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void trajectory<dimension, float_type>::depends()
{
    modules::depends<_Self, particle_type>::required();
}

template <int dimension, typename float_type>
trajectory<dimension, float_type>::trajectory(modules::factory& factory, po::options const& vm)
  // dependency injection
  : particle(modules::fetch<particle_type>(factory, vm))
  // allocate sample pointers
  , r(particle->ntype)
  , v(particle->ntype)
{
    for (size_t i = 0; i < particle->ntype; ++i) {
        r[i].reset(new sample_vector(particle->ntypes[i]));
        v[i].reset(new sample_vector(particle->ntypes[i]));
    }
}

#ifndef USE_HOST_SINGLE_PRECISION
template class trajectory<3, double>;
template class trajectory<2, double>;
#endif
template class trajectory<3, float>;
template class trajectory<2, float>;

}}} // namespace mdsim::samples::host

} // namespace halmd
