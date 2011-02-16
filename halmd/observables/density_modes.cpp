/*
 * Copyright © 2011  Felix Höfling
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

#include <boost/lexical_cast.hpp>
#include <string>

#include <halmd/observables/density_modes.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables
{

template <int dimension>
void density_modes<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("density_modes_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("observables")
            [
                class_<density_modes, shared_ptr<density_modes> >(class_name.c_str())
                    .def("acquire", &density_modes::acquire)
                    .property("wavenumbers", &density_modes::wavenumbers)
                    .property("ntype", &density_modes::ntype)
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__ ((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &density_modes<3>::luaopen
    ]
    [
        &density_modes<2>::luaopen
    ];
}

} // namespace

} // namespace observables

} // namespace halmd
