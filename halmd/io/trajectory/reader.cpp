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

#include <halmd/io/trajectory/reader.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace io { namespace trajectory
{

template <int dimension>
void reader<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("reader_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("io")
            [
                namespace_("trajectory")
                [
                    class_<reader, shared_ptr<reader> >(class_name.c_str())
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_io_trajectory_reader(lua_State* L)
{
    reader<3>::luaopen(L);
    reader<2>::luaopen(L);
    return 0;
}

}} // namespace io::trajectory

} // namespace halmd
