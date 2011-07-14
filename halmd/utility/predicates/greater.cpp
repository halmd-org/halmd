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

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include <halmd/utility/predicates/greater.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace predicates {

template <typename greater_type>
static typename greater_type::slot_function_type
wrap_evaluate(shared_ptr<greater_type const> greater)
{
    return bind(&greater_type::evaluate, greater);
}

template <typename value_type>
static shared_ptr<greater<value_type> >
wrap_greater(typename greater<value_type>::function_type const& func, value_type const& value)
{
    return make_shared<greater<value_type> >(func, value);
}

template <typename value_type>
void greater<value_type>::luaopen(lua_State* L, char const* class_name)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("predicates")
        [
            class_<greater>(class_name)
                .property("evaluate", &wrap_evaluate<greater>)
                .def("on_greater", &greater::on_greater)
                .scope
                [
                    class_<function_type>("function_type")
                        .def("__call", &function_type::operator())
                ]

          , def("greater", &wrap_greater<value_type>)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_utility_predicates_greater(lua_State* L)
{
    greater<double>::luaopen(L, "greater<double>");
    greater<float>::luaopen(L, "greater<float>");
    return 0;
}

} // namespace predicates
} // namespace halmd
