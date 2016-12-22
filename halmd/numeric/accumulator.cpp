/*
 * Copyright © 2016 Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <lua.hpp>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace numeric {

template <typename T>
void luaopen_accumulator(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("numeric")
        [
            class_<detail::accumulator<T>>()
                .def("sum", (T(*)(detail::accumulator<T> const&)) &detail::sum)
                .def("mean", (T const&(*)(detail::accumulator<T> const&)) &detail::mean)
                .def("error_of_mean", (T(*)(detail::accumulator<T> const&)) &detail::error_of_mean)
                .def("variance", (T(*)(detail::accumulator<T> const&)) &detail::variance)
                .def("count", (uint64_t const&(*)(detail::accumulator<T> const&)) &detail::count)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_numeric_accumulator(lua_State* L)
{
    luaopen_accumulator<double>(L);
    luaopen_accumulator<fixed_vector<double, 2>>(L);
    luaopen_accumulator<fixed_vector<double, 3>>(L);
    return 0;
}

} // namespace numeric
} // namespace halmd
