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

#include <halmd/mdsim/smoothers/localr4.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/demangle.hpp>

namespace halmd {
namespace mdsim {
namespace smoothers {

template <typename float_type>
void localr4<float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static std::string class_name("localr4_" + demangled_name<float_type>());
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("smoothers")
            [
                class_<localr4, boost::shared_ptr<localr4> >(class_name.c_str())
                    .def(constructor<float_type>())
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_smoothers_localr4(lua_State* L)
{
#if defined(HALMD_WITH_GPU) || defined(USE_HOST_SINGLE_PRECISION)
    localr4<float>::luaopen(L);
#endif
#ifndef USE_HOST_SINGLE_PRECISION
    localr4<double>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#if defined(HALMD_WITH_GPU) || defined(USE_HOST_SINGLE_PRECISION)
template class localr4<float>;
#endif
#ifndef USE_HOST_SINGLE_PRECISION
template class localr4<double>;
#endif

} // namespace smoothers
} // namespace mdsim
} // namespace halmd
