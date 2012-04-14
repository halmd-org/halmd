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
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
void force<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("force_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                class_<force, shared_ptr<_Base>, _Base>(class_name.c_str())
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_force(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    force<3, double>::luaopen(L);
    force<2, double>::luaopen(L);
#else
    force<3, float>::luaopen(L);
    force<2, float>::luaopen(L);
#endif
    return 0;
}

} // namespace mdsim
} // namespace host
} // namespace halmd