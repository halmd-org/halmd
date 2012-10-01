/*
 * Copyright © 2011  Felix Höfling
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
#include <boost/bind.hpp>

#include <halmd/observables/fields/density.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;

namespace halmd {
namespace observables {
namespace fields {

template <typename density_type>
static typename density_type::slot_function_type
sample_wrapper(shared_ptr<density_type> density)
{
    return bind(&density_type::sample, density);
}

template <typename density_type>
static function<typename density_type::result_type const& ()>
wrap_value(shared_ptr<density_type> density)
{
    return bind(&density_type::value, density);
}

template <int dimension>
void density<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("fields")
            [
                class_<density, shared_ptr<density> >("density_")
                    .property("sample", &sample_wrapper<density>)
                    .property("value", &wrap_value<density>)
                    .def("on_sample", &density::on_sample)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_fields_density(lua_State* L)
{
    density<2>::luaopen(L);
    density<3>::luaopen(L);
    return 0;
}

} // namespace fields
} // namespace observables
} // namespace halmd
