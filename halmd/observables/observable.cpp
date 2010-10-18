/*
 * Copyright © 2010  Felix Höfling
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

#include <halmd/observables/observable.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables
{

template <int dimension>
void observable<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    string class_name("observable_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            class_<observable, shared_ptr<observable> >(class_name.c_str())
                .def("register_observables", &observable::register_observables)
                .def("sample", &observable::sample)
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &observable<3>::luaopen
    ]
    [
        &observable<2>::luaopen
    ];
}

} // namespace observables

} // namespace halmd
