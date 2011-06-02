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

#include <halmd/mdsim/position.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim
{

template <typename position>
typename signal<void ()>::slot_function_type
wrap_set(shared_ptr<position> self)
{
    return bind(&position::set, self);
}

template <int dimension>
void position<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("position_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<position, shared_ptr<position> >(class_name.c_str())
                .property("set", &wrap_set<position>)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_position(lua_State* L)
{
    position<3>::luaopen(L);
    position<2>::luaopen(L);
    return 0;
}

template class position<3>;
template class position<2>;

} // namespace mdsim

} // namespace halmd
