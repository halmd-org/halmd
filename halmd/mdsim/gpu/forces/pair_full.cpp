/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/forces/pair_full.hpp>
// #include <halmd/mdsim/gpu/forces/smooth.hpp>
#include <halmd/mdsim/gpu/potentials/lennard_jones.hpp>
#include <halmd/mdsim/gpu/potentials/lennard_jones_simple.hpp>
#include <halmd/mdsim/gpu/potentials/modified_lennard_jones.hpp>
#include <halmd/mdsim/gpu/potentials/morse.hpp>
#include <halmd/mdsim/gpu/potentials/power_law.hpp>
#include <halmd/mdsim/gpu/potentials/power_law_with_core.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

using namespace potentials;

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_forces_pair_full(lua_State* L)
{
    pair_full<3, float, lennard_jones<float> >::luaopen(L);
    pair_full<2, float, lennard_jones<float> >::luaopen(L);

    pair_full<3, float, lennard_jones_simple<float> >::luaopen(L);
    pair_full<2, float, lennard_jones_simple<float> >::luaopen(L);

    pair_full<3, float, modified_lennard_jones<float> >::luaopen(L);
    pair_full<2, float, modified_lennard_jones<float> >::luaopen(L);

    pair_full<3, float, morse<float> >::luaopen(L);
    pair_full<2, float, morse<float> >::luaopen(L);

    pair_full<3, float, power_law<float> >::luaopen(L);
    pair_full<2, float, power_law<float> >::luaopen(L);

    pair_full<3, float, power_law_with_core<float> >::luaopen(L);
    pair_full<2, float, power_law_with_core<float> >::luaopen(L);

    return 0;
}

// explicit instantiation of force modules
template class pair_full<3, float, lennard_jones<float> >;
template class pair_full<2, float, lennard_jones<float> >;

template class pair_full<3, float, lennard_jones_simple<float> >;
template class pair_full<2, float, lennard_jones_simple<float> >;

template class pair_full<3, float, modified_lennard_jones<float> >;
template class pair_full<2, float, modified_lennard_jones<float> >;

template class pair_full<3, float, morse<float> >;
template class pair_full<2, float, morse<float> >;

template class pair_full<3, float, power_law<float> >;
template class pair_full<2, float, power_law<float> >;

template class pair_full<3, float, power_law_with_core<float> >;
template class pair_full<2, float, power_law_with_core<float> >;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
