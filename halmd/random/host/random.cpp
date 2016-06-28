/*
 * Copyright © 2013-2014 Felix Höfling
 * Copyright © 2008-2011 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <memory>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace random {
namespace host {

std::shared_ptr<logger> const logger_ = std::make_shared<logger>("random (host)");

random::random(unsigned int seed)
{
    LOG("random number generator type: " << rng_name());
    random::seed(seed);
}

void random::seed(unsigned int seed)
{
    LOG("set RNG seed: " << seed);
    rng_.seed(seed);
}

/**
 * shuffle sequence of abstract Lua type (number, table, userdata, …)
 * and returned shuffled sequence
 */
static std::vector<luaponte::object>
wrap_shuffle(std::shared_ptr<random> self, std::vector<luaponte::object> const& v)
{
    std::vector<luaponte::object> result(v);    // FIXME deep copy needed due to luaponte::out_value bug
    self->shuffle(result.begin(), result.end());
    return result;
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
                class_<random, std::shared_ptr<random>>(rng_name())
                    .def(constructor<>())
                    .def(constructor<unsigned int>())
                    .def("seed", &random::seed)
                    .def("shuffle", &wrap_shuffle)
//                    .def("shuffle", &wrap_shuffle, out_value(_2)) FIXME does not compile
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_random_host_random(lua_State* L)
{
    random::luaopen(L);
    return 0;
}

} // namespace host
} // namespace random
} // namespace halmd
