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
#include <boost/iterator/counting_iterator.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/velocities/trajectory.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace velocities
{

template <int dimension, typename float_type>
trajectory<dimension, float_type>::trajectory(
    shared_ptr<particle_type> particle
  , shared_ptr<sample_type> sample
)
  : _Base(particle)
  // dependency injection
  , particle(particle)
  , sample(sample)
{
}

/**
 * set particle velocities
 */
template <int dimension, typename float_type>
void trajectory<dimension, float_type>::set()
{
    for (size_t j = 0, i = 0; j < particle->ntype; i += particle->ntypes[j], ++j) {
        copy(sample->v[j]->begin(), sample->v[j]->end(), &particle->v[i]);
    }

    LOG("set particle velocities from trajectory sample");
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class trajectory<3, double>;
template class trajectory<2, double>;
#else
template class trajectory<3, float>;
template class trajectory<2, float>;
#endif

}}} // namespace mdsim::host::velocities

} // namespace halmd
