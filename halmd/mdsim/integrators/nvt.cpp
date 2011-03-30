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
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

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
    module(L)
    [
        namespace_("libhalmd")
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
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &nvt<3>::luaopen
    ]
    [
        &nvt<2>::luaopen
    ];
}

} // namespace
// explicit instantiation
template class nvt<3>;
template class nvt<2>;

}} // namespace mdsim::integrators

} // namespace halmd
