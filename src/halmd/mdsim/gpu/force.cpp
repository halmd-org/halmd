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

#include <halmd/mdsim/gpu/force.hpp>

using namespace boost;

namespace halmd
{
namespace mdsim { namespace gpu
{

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void force<dimension, float_type>::depends()
{
    modules::depends<_Self, particle_type>::required();
    modules::depends<_Self, box_type>::required();
//     modules::depends<_Self, smooth_type>::optional();
}

/**
 * Initialize module dependencies
 */
template <int dimension, typename float_type>
force<dimension, float_type>::force(modules::factory& factory, po::variables_map const& vm)
  : _Base(factory, vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(factory, vm))
  , box(modules::fetch<box_type>(factory, vm))
//   , smooth(modules::fetch<smooth_type>(factory, vm))
  // allocate result variables
  , g_en_pot_(particle->dim.threads())
  , g_stress_pot_(particle->dim.threads())
{
}

// explicit instantiation
template class force<3, float>;
template class force<2, float>;

}} // namespace mdsim::gpu

} // namespace halmd
