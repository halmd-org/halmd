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

#include <boost/bind.hpp>

#include <halmd/io/trajectory/writer.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace io {
namespace trajectory {

template <typename writer_type>
typename signal<void (uint64_t)>::slot_function_type
append_wrapper(shared_ptr<writer_type> writer)
{
    return bind(&writer_type::append, writer, _1);
}

template <int dimension>
void writer<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("writer_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            namespace_("trajectory")
            [
                class_<writer, shared_ptr<writer> >(class_name.c_str())
                    .property("append", &append_wrapper<writer>)
                    .def("on_append", &writer::on_append)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_trajectory_writer(lua_State* L)
{
    writer<3>::luaopen(L);
    writer<2>::luaopen(L);
    return 0;
}

} // namespace io
} // namespace trajectory
} // namespace halmd
