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

#include <boost/lexical_cast.hpp>
#include <string>

#include <halmd/observables/gpu/samples/phase_space.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd {
namespace observables {
namespace gpu {
namespace samples {

template <int dimension, typename float_type>
char const* phase_space<dimension, float_type>::class_name()
{
    static string class_name("phase_space_" + lexical_cast<string>(dimension) + "_");
    return class_name.c_str();
}

template <int dimension, typename float_type>
static char const* class_name_wrapper(phase_space<dimension, float_type> const&)
{
    return phase_space<dimension, float_type>::class_name();
}

template <int dimension, typename float_type>
static int wrap_dimension(phase_space<dimension, float_type> const&)
{
    return dimension;
}

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                namespace_("samples")
                [
                    class_<phase_space, shared_ptr<phase_space> >(class_name())
                        .def(constructor<vector<unsigned int> >())
                        .property("class_name", &class_name_wrapper<dimension, float_type>)
                        .property("dimension", &wrap_dimension<dimension, float_type>)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_samples_phase_space(lua_State* L)
{
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
    observables::samples::blocking_scheme<phase_space<3, float> >::luaopen(L, "gpu");
    observables::samples::blocking_scheme<phase_space<2, float> >::luaopen(L, "gpu");
    return 0;
}

// explicit instantiation
template class phase_space<3, float>;
template class phase_space<2, float>;

} // namespace samples
} // namespace gpu

namespace samples
{

// explicit instantiation
template class blocking_scheme<gpu::samples::phase_space<3, float> >;
template class blocking_scheme<gpu::samples::phase_space<2, float> >;

} // namespace samples
} // namespace observables
} // namespace halmd
