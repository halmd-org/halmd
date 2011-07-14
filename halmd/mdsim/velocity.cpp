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

#include <halmd/mdsim/velocity.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

template <typename velocity>
typename signal<void ()>::slot_function_type
wrap_set(shared_ptr<velocity> self)
{
    return bind(&velocity::set, self);
}

template <int dimension>
void velocity<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("velocity_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<velocity, shared_ptr<velocity> >(class_name.c_str())
                .property("set", &wrap_set<velocity>)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_velocity(lua_State* L)
{
    velocity<3>::luaopen(L);
    velocity<2>::luaopen(L);
    return 0;
}

template class velocity<3>;
template class velocity<2>;

} // namespace mdsim
} // namespace halmd
