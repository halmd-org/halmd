/*
 * Copyright © 2011-2012 Peter Colberg
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

#include <halmd/observables/samples/blocking_scheme.hpp>

namespace halmd {
namespace observables {
namespace samples {

/**
 * This function transforms a Lua function into a std::function, for
 * use with the constructor of blocking_scheme<luaponte::object>.
 */
static std::function<std::shared_ptr<luaponte::object const> ()>
blocking_scheme_adaptor(luaponte::object const& function)
{
    using namespace luaponte;
    return [=]() {
        return std::make_shared<object>(call_function<object>(function));
    };
}

HALMD_LUA_API int luaopen_libhalmd_observables_samples_blocking_scheme(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("samples")
            [
                class_<blocking_scheme_base>()

              , def("blocking_scheme_adaptor", &blocking_scheme_adaptor)
            ]
        ]
    ];
    blocking_scheme<luaponte::object>::luaopen(L);
    return 0;
}

// explicit instantiation
template class blocking_scheme<luaponte::object>;

} // namespace samples
} // namespace observables
} // namespace halmd
