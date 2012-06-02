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

#include <halmd/observables/samples/blocking_scheme.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace samples {

HALMD_LUA_API int luaopen_libhalmd_observables_samples_blocking_scheme(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("samples")
            [
                class_<blocking_scheme_base, std::shared_ptr<blocking_scheme_base> >("blocking_scheme_base")
            ]
        ]
    ];
    return 0;
}

} // namespace samples
} // namespace observables
} // namespace halmd
