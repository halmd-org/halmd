/*
 * Copyright © 2010-2011  Peter Colberg
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
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace utility
{

profiler::profiler(writers_type writer, string const& group)
  : writer_(writer)
  , group_(group)
{
    LOG_ONCE("timer resolution: " << 1.E9 * timer::elapsed_min() << " ns");
}

void profiler::register_runtime(accumulator_type const& runtime, string const& tag, std::string const& desc)
{
    writers_type::const_iterator i, ie = writer_.end();
    for (i = writer_.begin(); i != ie; ++i) {
        (*i)->register_accumulator(group_, runtime, tag, desc);
    }
}

void profiler::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("utility")
            [
                class_<profiler, shared_ptr<profiler> >("profiler")
                    .def(constructor<writers_type, string>())
            ]
        ]
    ];
}

HALMD_INIT( register_luaopen )
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &profiler::luaopen
    ];
}

} // namespace utility

} // namespace halmd
