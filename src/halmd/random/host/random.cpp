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

template <typename Module>
static void register_lua(char const* class_name)
{
    typedef typename Module::_Base _Base;

    using namespace luabind;
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        namespace_("halmd_wrapper")
        [
            namespace_("host")
            [
                namespace_("random")
                [
                    class_<Module, shared_ptr<_Base>, bases<_Base> >(class_name)
                        .def(constructor<unsigned int>())
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    register_lua<random>("gfsr4");
}

}} // namespace random::host

} // namespace halmd
