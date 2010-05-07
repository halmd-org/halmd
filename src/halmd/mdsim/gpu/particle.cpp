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
#include <exception>
#include <numeric>

#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/util/logger.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu
{

/**
 * Resolve module dependencies
 */
template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::resolve(po::options const& vm)
{
    if (vm["backend"].as<string>() != "gpu_neighbour") {
        throw module_error("inept module " + module<particle>::name());
    }
    module<device_type>::required(vm);
}

template <unsigned int dimension, typename float_type>
particle<dimension, float_type>::particle(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , device(module<device_type>::fetch(vm))
  // allocate global device memory
  , g_r(nbox)
  , g_image(nbox)
  , g_v(nbox)
  , g_f(nbox)
  , g_neighbor(nbox)
  // allocate page-locked host memory
  , h_r(nbox)
  , h_image(nbox)
  , h_v(nbox)
{
}

// explicit instantiation
template class particle<3, float>;
template class particle<2, float>;

}} // namespace mdsim::gpu

template class module<mdsim::gpu::particle<3, float> >;
template class module<mdsim::gpu::particle<2, float> >;

} // namespace halmd
