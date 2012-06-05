/*
 * Copyright © 2012 Peter Colberg
 * Copyright © 2012 Felix Höfling
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

#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {

void particle_group::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L)
    [
        class_<particle_group>()
            .property("size", &particle_group::size)
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_particle_group(lua_State* L)
{
    particle_group::luaopen(L);
    return 0;
}

} // namespace host
} // namespace mdsim
} // namespace halmd
