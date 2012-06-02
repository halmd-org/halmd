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

#include <memory>

#include <halmd/utility/predicates/greater.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace predicates {

template <typename greater_type>
static typename greater_type::slot_function_type
wrap_evaluate(std::shared_ptr<greater_type const> self)
{
    return [=]() {
        self->evaluate();
    };
}

template <typename value_type>
void greater<value_type>::luaopen(lua_State* L, char const* class_name)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("predicates")
        [
            class_<greater>(class_name)
                .property("evaluate", &wrap_evaluate<greater>)
                .def("on_greater", &greater::on_greater)

          , def("greater", &std::make_shared<greater, function_type, value_type>)
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
