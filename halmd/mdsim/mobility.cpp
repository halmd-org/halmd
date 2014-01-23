/*
 * Copyright © 2011 Michael Kopp
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

#include <halmd/mdsim/mobility.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

template <typename mobility_type>
typename signal<void ()>::slot_function_type
wrap_compute(shared_ptr<mobility_type> mobility)
{
    return bind(&mobility_type::compute, mobility);
}

template <typename mobility_type>
typename signal<void ()>::slot_function_type
wrap_compute_velocity(shared_ptr<mobility_type> mobility)
{
    return bind(&mobility_type::compute_velocity, mobility);
}

template <int dimension>
void mobility<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("mobility_" + lexical_cast<string>(dimension) + "_");

    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<mobility, shared_ptr<mobility> >(class_name.c_str())
                .def("compute", &mobility::compute)
                .def("compute_velocity", &mobility::compute_velocity)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_mobility(lua_State* L)
{
    mobility<3>::luaopen(L);
    mobility<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class mobility<3>;
template class mobility<2>;

} // namespace mdsim
} // namespace halmd
