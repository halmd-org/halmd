/*
 * Copyright © 2012  Nicolas Höft
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

#include <halmd/mdsim/smoothers/nosmooth.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace smoothers {

void nosmooth::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("smoothers")
            [
                class_<nosmooth, boost::shared_ptr<nosmooth> >("nosmooth")
                    .def(constructor<>())
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_smoothers_nosmooth(lua_State* L)
{
    nosmooth::luaopen(L);
    return 0;
}

} // namespace smoothers
} // namespace mdsim
} // namespace halmd
