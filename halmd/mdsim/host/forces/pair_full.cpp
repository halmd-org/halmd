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

#include <halmd/mdsim/host/forces/pair_full.hpp>
#include <halmd/mdsim/host/forces/smooth.hpp>
#include <halmd/mdsim/host/potentials/lennard_jones.hpp>
#include <halmd/mdsim/host/potentials/modified_lennard_jones.hpp>
#include <halmd/mdsim/host/potentials/morse.hpp>
#include <halmd/mdsim/host/potentials/power_law.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace forces {

using namespace potentials;

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_forces_pair_full(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    typedef double float_type;
#else
    typedef float float_type;
#endif

    pair_full<3, float_type, lennard_jones<float_type> >::luaopen(L);
    pair_full<2, float_type, lennard_jones<float_type> >::luaopen(L);

    pair_full<3, float_type, modified_lennard_jones<float_type> >::luaopen(L);
    pair_full<2, float_type, modified_lennard_jones<float_type> >::luaopen(L);

    pair_full<3, float_type, morse<float_type> >::luaopen(L);
    pair_full<2, float_type, morse<float_type> >::luaopen(L);

    pair_full<3, float_type, power_law<float_type> >::luaopen(L);
    pair_full<2, float_type, power_law<float_type> >::luaopen(L);
    return 0;
}

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
typedef double float_type;
#else
typedef float float_type;
#endif

template class pair_full<3, float_type, lennard_jones<float_type> >;
template class pair_full<2, float_type, lennard_jones<float_type> >;

template class pair_full<3, float_type, modified_lennard_jones<float_type> >;
template class pair_full<2, float_type, modified_lennard_jones<float_type> >;

template class pair_full<3, float_type, morse<float_type> >;
template class pair_full<2, float_type, morse<float_type> >;

template class pair_full<3, float_type, power_law<float_type> >;
template class pair_full<2, float_type, power_law<float_type> >;

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
