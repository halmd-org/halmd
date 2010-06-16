/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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
#include <halmd/mdsim/gpu/integrator.hpp>
#include <halmd/utility/module.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu
{

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void integrator<dimension, float_type>::options(po::options_description& desc)
{
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void integrator<dimension, float_type>::depends()
{
    modules::required<_Self, particle_type>();
    modules::required<_Self, box_type>();
    modules::required<_Self, device_type>();
}

template <int dimension, typename float_type>
integrator<dimension, float_type>::integrator(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(vm))
  , box(modules::fetch<box_type>(vm))
  , device(modules::fetch<device_type>(vm))
{
}

// explicit instantiation
template class integrator<3, float>;
template class integrator<2, float>;

}} // namespace mdsim::gpu

} // namespace halmd
