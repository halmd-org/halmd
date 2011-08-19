/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/force.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

template <typename force_type>
typename signal<void ()>::slot_function_type
wrap_compute(shared_ptr<force_type> force)
{
    return bind(&force_type::compute, force);
}

template <typename force_type>
typename signal<void ()>::slot_function_type
wrap_aux_enable(shared_ptr<force_type> force)
{
    return bind(&force_type::aux_enable, force);
}

template <typename force_type>
typename signal<void ()>::slot_function_type
wrap_aux_disable(shared_ptr<force_type> force)
{
    return bind(&force_type::aux_disable, force);
}

template <int dimension>
void force<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("force_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<force, shared_ptr<force> >(class_name.c_str())
                .property("compute", &wrap_compute<force>)
                .property("aux_enable", &wrap_aux_enable<force>)
                .property("aux_disable", &wrap_aux_disable<force>)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_force(lua_State* L)
{
    force<3>::luaopen(L);
    force<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class force<3>;
template class force<2>;

} // namespace mdsim
} // namespace halmd
