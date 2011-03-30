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

#include <halmd/io/profiling/writer.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace profiling
{

template <typename writer_type>
typename signal<void (double)>::slot_function_type
write_wrapper(shared_ptr<writer_type> writer)
{
    return bind(&writer_type::write, writer);
}

void writer::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("io")
            [
                namespace_("profiling")
                [
                    class_<writer, shared_ptr<writer> >("writer")
                        .property("write", &write_wrapper<writer>)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_profiling_writer(lua_State* L)
{
    writer::luaopen(L);
    return 0;
}

}} // namespace io::profiling

} // namespace halmd

