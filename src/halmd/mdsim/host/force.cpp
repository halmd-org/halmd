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

#include <halmd/mdsim/host/force.hpp>

using namespace boost;

namespace halmd
{
namespace mdsim { namespace host
{

/**
 * Assemble module options
 */
template <int dimension, typename float_type>
void force<dimension, float_type>::options(po::options_description& desc)
{
}

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void force<dimension, float_type>::resolve(po::options const& vm)
{
    module<particle_type>::required(vm);
    module<box_type>::required(vm);
    module<neighbour_type>::required(vm);
    module<thermodynamics_type>::required(vm);
    module<smooth_type>::optional(vm);
}

/**
 * Initialize module dependencies
 */
template <int dimension, typename float_type>
force<dimension, float_type>::force(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  , box(module<box_type>::fetch(vm))
  , neighbour(module<neighbour_type>::fetch(vm))
  , thermodynamics(module<thermodynamics_type>::fetch(vm))
  , smooth(module<smooth_type>::fetch(vm))
{
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class force<3, double>;
template class force<2, double>;
#else
template class force<3, float>;
template class force<2, float>;
#endif

}} // namespace mdsim::host

} // namespace halmd
