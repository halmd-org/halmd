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
#include <halmd/mdsim/integrator.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace mdsim {

template <typename integrator_type>
typename signal<void ()>::slot_function_type
wrap_integrate(shared_ptr<integrator_type> integrator)
{
    return bind(&integrator_type::integrate, integrator);
}

template <typename integrator_type>
typename signal<void ()>::slot_function_type
wrap_finalize(shared_ptr<integrator_type> integrator)
{
    return bind(&integrator_type::finalize, integrator);
}

template <int dimension>
void integrator<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("integrator_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<integrator, shared_ptr<integrator> >(class_name.c_str())
                .property("timestep", static_cast<double (integrator::*)() const>(&integrator::timestep), static_cast<void (integrator::*)(double)>(&integrator::timestep))
                .property("integrate", &wrap_integrate<integrator>)
                .property("finalize", &wrap_finalize<integrator>)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_integrator(lua_State* L)
{
    integrator<3>::luaopen(L);
    integrator<2>::luaopen(L);
    return 0;
}

// explicit instantiation
template class integrator<3>;
template class integrator<2>;

} // namespace mdsim
} // namespace halmd
