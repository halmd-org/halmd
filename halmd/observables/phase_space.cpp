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

#include <halmd/observables/phase_space.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/signal.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables
{

template <typename phase_space_type>
typename signal<void (double)>::slot_function_type
acquire_wrapper(shared_ptr<phase_space_type> phase_space)
{
    return bind(&phase_space_type::acquire, phase_space, _1);
}

template <int dimension>
void phase_space<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("phase_space_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            class_<phase_space, shared_ptr<phase_space> >(class_name.c_str())
                .property("acquire", &acquire_wrapper<phase_space>)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_phase_space(lua_State* L)
{
    phase_space<3>::luaopen(L);
    phase_space<2>::luaopen(L);
    return 0;
}

} // namespace observables

} // namespace halmd
