/*
 * Copyright © 2010  Peter Colberg
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

#include <halmd/random/random.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/read_integer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace random
{

void random::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("random")
            [
                class_<random, shared_ptr<random> >("random")
                    .scope
                    [
                        def("read_integer", &read_integer<unsigned int>)
                    ]
            ]
        ]
    ];
}

namespace // limit symbols to translation unit
{

__attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &random::luaopen
    ];
}

} // namespace

} // namespace random

} // namespace halmd
