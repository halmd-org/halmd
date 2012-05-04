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

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/sort.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

template <typename sort_type>
typename signal<void ()>::slot_function_type
wrap_order(boost::shared_ptr<sort_type> sort)
{
    return bind(&sort_type::order, sort);
}

template <int dimension>
void sort<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("sort_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<sort, boost::shared_ptr<sort> >(class_name.c_str())
                .property("order", &wrap_order<sort>)
                .def("on_order", &sort::on_order)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_sort(lua_State* L)
{
    sort<3>::luaopen(L);
    sort<2>::luaopen(L);
    return 0;
}

} // namespace mdsim
} // namespace halmd
