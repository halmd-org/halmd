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
#include <halmd/random/host/random.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/read_integer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace random { namespace host
{

random::random(unsigned int seed)
{
    LOG("random number generator seed: " << seed);
    rng_.seed(seed);
}

void random::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("host")
            [
                namespace_("random")
                [
                    class_<random, shared_ptr<_Base>, bases<_Base> >("gfsr4")
                        .def(constructor<unsigned int>())
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &random::luaopen
    ];
}

}} // namespace random::host

} // namespace halmd
