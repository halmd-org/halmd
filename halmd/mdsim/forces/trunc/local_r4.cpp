/*
 * Copyright © 2008-2010 Peter Colberg and Felix Höfling
 * Copyright © 2012 Nicolas Höft
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

#include <string>

#include <halmd/mdsim/forces/trunc/local_r4.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/demangle.hpp>

namespace halmd {
namespace mdsim {
namespace forces {
namespace trunc {

template <typename float_type>
void local_r4<float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static std::string class_name("local_r4_" + demangled_name<float_type>());
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("forces")
            [
                namespace_("trunc")
                [
                    class_<local_r4, boost::shared_ptr<local_r4> >(class_name.c_str())
                        .def(constructor<float_type>())
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_forces_trunc_local_r4(lua_State* L)
{
    local_r4<float>::luaopen(L);
    local_r4<double>::luaopen(L);
    return 0;
}

// explicit instantiation
template class local_r4<float>;
template class local_r4<double>;

} // namespace trunc
} // namespace forces
} // namespace mdsim
} // namespace halmd
