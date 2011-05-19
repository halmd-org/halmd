/*
 * Copyright © 2011  Peter Colberg
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

#include <boost/shared_ptr.hpp>

#include <halmd/observables/dynamics/blocking_scheme.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace dynamics
{

blocking_scheme::blocking_scheme(shared_ptr<block_sample_type> block_sample)
  : block_sample_(block_sample)
{
}

HALMD_LUA_API int luaopen_libhalmd_observables_dynamics_blocking_scheme(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<blocking_scheme, shared_ptr<blocking_scheme> >("blocking_scheme")
                    .def(constructor<shared_ptr<blocking_scheme::block_sample_type> >())
            ]
        ]
    ];
    return 0;
}

}} // namespace observables::dynamics

} // namespace halmd
