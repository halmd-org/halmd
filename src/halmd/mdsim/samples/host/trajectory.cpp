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

#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/sample/trajectory.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
#endif
#include <halmd/mdsim/host/sample/trajectory.hpp>
#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace mdsim { namespace samples { namespace host
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void trajectory<dimension, float_type>::resolve(po::options const& vm)
{
    module<particle_type>::required(vm);
}

template <int dimension, typename float_type>
trajectory<dimension, float_type>::trajectory(po::options const& vm)
  // dependency injection
  : particle(module<particle_type>::fetch(vm))
  // allocate sample pointers
  , r(particle->ntype)
  , v(particle->ntype)
{}

#ifndef USE_HOST_SINGLE_PRECISION
template class trajectory<3, double>;
template class trajectory<2, double>;
#endif
template class trajectory<3, float>;
template class trajectory<2, float>;

}}} // namespace mdsim::samples::host

} // namespace halmd
