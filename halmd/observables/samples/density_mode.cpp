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

#include <boost/lexical_cast.hpp>
#include <string>

#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/observables/samples/density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace samples {

template <int dimension>
char const* density_mode<dimension>::class_name()
{
    static string class_name("density_mode_" + lexical_cast<string>(dimension) + "_");
    return class_name.c_str();
}

template <int dimension>
static char const* class_name_wrapper(density_mode<dimension> const&)
{
    return density_mode<dimension>::class_name();
}

template <int dimension>
void density_mode<dimension>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("samples")
            [
                class_<density_mode, shared_ptr<density_mode> >(class_name())
                    .def(constructor<unsigned int, unsigned int>())
                    .property("class_name", &class_name_wrapper<dimension>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_samples_density_mode(lua_State* L)
{
    density_mode<3>::luaopen(L);
    density_mode<2>::luaopen(L);
    observables::samples::blocking_scheme<density_mode<3> >::luaopen(L);
    observables::samples::blocking_scheme<density_mode<2> >::luaopen(L);
    return 0;
}

template class density_mode<3>;
template class density_mode<2>;
template class blocking_scheme<density_mode<3> >;
template class blocking_scheme<density_mode<2> >;

} // namespace samples
} // namespace observables
} // namespace halmd
