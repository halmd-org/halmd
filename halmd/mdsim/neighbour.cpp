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

#include <boost/bind.hpp>

#include <halmd/mdsim/neighbour.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

template <typename neighbour_type>
typename signal<void ()>::slot_function_type
wrap_update(shared_ptr<neighbour_type> neighbour)
{
    return bind(&neighbour_type::update, neighbour);
}

void neighbour::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<neighbour>("neighbour")
                .property("update", &wrap_update<neighbour>)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_neighbour(lua_State* L)
{
    neighbour::luaopen(L);
    return 0;
}

} // namespace mdsim
} // namespace halmd
