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

#include <boost/nondet_random.hpp> // boost::random_device

#include <halmd/io/logger.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace random {
namespace host {

random::random(unsigned int seed)
{
    LOG("random number generator seed: " << seed);
    rng_.seed(seed);
}

//! Get seed from non-deterministic random number generator.
// boost::random_device reads from /dev/urandom on GNU/Linux,
// and the default cryptographic service provider on Windows.
unsigned int random::defaults::seed() {
    return boost::random_device()();
}

void random::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("host")
        [
            namespace_("random")
            [
                class_<random, shared_ptr<random> >("gfsr4")
                    .def(constructor<unsigned int>())
                    .scope
                    [
                        class_<defaults>("defaults")
                            .scope
                            [
                                def("seed", &defaults::seed)
                            ]
                    ]
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
