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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/neighbour.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

template <int dimension>
void neighbour<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("neighbour_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                class_<neighbour, shared_ptr<neighbour> >(class_name.c_str())
                    .def("check", &neighbour::check)
                    .def("update", &neighbour::update)
            ]
        ]
    ];
}

HALMD_INIT( register_luaopen )
{
    lua_wrapper::register_(0) //< distance to base class
    [
        &neighbour<3>::luaopen
    ]
    [
        &neighbour<2>::luaopen
    ];
}

template class neighbour<3>;
template class neighbour<2>;

} // namespace mdsim

} // namespace halmd
