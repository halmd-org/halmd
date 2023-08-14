/*
 * Copyright © 2014-2015 Nicolas Höft
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/mdsim/gpu/particle_groups/region_kernel.cu>
#include <test/unit/mdsim/geometries/simple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_groups{

template class region_wrapper<2, simple_geometry<2, float> >;
template class region_wrapper<3, simple_geometry<3, float> >;

} // namespace particle_groups
} // namespace gpu
} // namespace mdsim
} // namespace halmd
