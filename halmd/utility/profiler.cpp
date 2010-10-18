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

#include <halmd/io/logger.hpp>
#include <halmd/io/profile/writer.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace utility
{

profiler::profiler(vector<shared_ptr<profile_writer_type> > profile_writers)
  // dependency injection
  : profile_writers(profile_writers)
{
    LOG("timer resolution: " << 1.E9 * timer::elapsed_min() << " ns");
}

void profiler::register_accumulator(
    std::type_info const& tag
  , accumulator_type const& acc
  , std::string const& desc
) const
{
    for_each(
        profile_writers.begin()
      , profile_writers.end()
      , bind(
            &profile_writer_type::register_accumulator
          , _1
          , tokenized_name(tag) // namespace and class template tokens
          , cref(acc)
          , cref(desc)
        )
    );
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
                    .def(constructor<vector<shared_ptr<profile_writer_type> > >())
                    .def_readonly("profile_writers", &profiler::profile_writers)
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &profiler::luaopen
    ];
}

} // namespace utility

} // namespace halmd
