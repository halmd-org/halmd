/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <memory>

#include <halmd/random/host/random.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace random {
namespace host {

random::random(unsigned int seed)
{
    random::seed(seed);
}

void random::seed(unsigned int seed)
{
    rng_.seed(seed);
}

void random::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("random")
        [
            namespace_("host")
            [
                class_<random, std::shared_ptr<random> >("gfsr4")
                    .def(constructor<>())
                    .def("seed", &random::seed)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_random_host_random(lua_State* L)
{
    random::luaopen(L);
    return 0;
}

} // namespace random
} // namespace host
} // namespace halmd
