/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <limits>

#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace host { namespace samples
{

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("phase_space_" + lexical_cast<string>(dimension) + "_" + demangled_name<float_type>() + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                namespace_("samples")
                [
                    class_<phase_space, shared_ptr<phase_space> >(class_name.c_str())
                        .def(constructor<vector<unsigned int> >())
                ]
            ]
        ]
    ];
    observables::samples::blocking_scheme<phase_space>::luaopen(L, "host", class_name.c_str());
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_samples_phase_space(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    phase_space<3, double>::luaopen(L);
    phase_space<2, double>::luaopen(L);
#endif
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class phase_space<3, double>;
template class phase_space<2, double>;
#endif
template class phase_space<3, float>;
template class phase_space<2, float>;

}}} // namespace observables::host::samples

namespace observables { namespace samples
{

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class blocking_scheme<observables::host::samples::phase_space<3, double> >;
template class blocking_scheme<observables::host::samples::phase_space<2, double> >;
#endif
template class blocking_scheme<observables::host::samples::phase_space<3, float> >;
template class blocking_scheme<observables::host::samples::phase_space<2, float> >;

}} // namespace observables::samples

} // namespace halmd
