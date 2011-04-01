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

#include <halmd/mdsim/integrators/nvt.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace integrators
{

template <int dimension>
void nvt<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("nvt_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("integrators")
            [
                class_<nvt, shared_ptr<_Base>, bases<_Base> >(class_name.c_str())
                    .property(
                        "temperature"
                      , (double (nvt::*)() const) &nvt::temperature
                    )
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_integrators_nvt(lua_State* L)
{
    nvt<3>::luaopen(L);
    nvt<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class nvt<3>;
template class nvt<2>;

}} // namespace mdsim::integrators

} // namespace halmd
