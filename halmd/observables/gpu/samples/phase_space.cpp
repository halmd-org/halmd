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

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
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
static int wrap_dimension(phase_space<dimension, float_type> const&)
{
    return dimension;
}

#ifndef NDEBUG
template <int dimension, typename float_type>
vector<typename phase_space<dimension, float_type>::vector_type>
phase_space<dimension, float_type>::position(unsigned int type) const
{
    if (!(type < r.size())) {
        throw invalid_argument("particle type");
    }
    cuda::vector<gpu_vector_type> const& g_r = *r[type];
    cuda::host::vector<gpu_vector_type> h_r(g_r.size());
    cuda::copy(g_r, h_r);
    vector<vector_type> position(g_r.size());
    copy(h_r.begin(), h_r.end(), position.begin());
    return position;
}

template <int dimension, typename float_type>
vector<typename phase_space<dimension, float_type>::vector_type>
phase_space<dimension, float_type>::velocity(unsigned int type) const
{
    if (!(type < v.size())) {
        throw invalid_argument("particle type");
    }
    cuda::vector<gpu_vector_type> const& g_v = *v[type];
    cuda::host::vector<gpu_vector_type> h_v(g_v.size());
    cuda::copy(g_v, h_v);
    vector<vector_type> velocity(g_v.size());
    copy(h_v.begin(), h_v.end(), velocity.begin());
    return velocity;
}

template <typename phase_space_type>
static function<vector<typename phase_space_type::vector_type> ()>
wrap_position(boost::shared_ptr<phase_space_type> self, unsigned int type)
{
    return bind(&phase_space_type::position, self, type);
}

template <typename phase_space_type>
static function<vector<typename phase_space_type::vector_type> ()>
wrap_velocity(boost::shared_ptr<phase_space_type> self, unsigned int type)
{
    return bind(&phase_space_type::velocity, self, type);
}
#endif /* ! NDEBUG */

template <int dimension, typename float_type>
void phase_space<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string const class_name("phase_space_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                namespace_("samples")
                [
                    class_<phase_space, boost::shared_ptr<phase_space> >(class_name.c_str())
                        .def(constructor<vector<unsigned int> >())
                        .property("dimension", &wrap_dimension<dimension, float_type>)
#ifndef NDEBUG
                        .def("position", &wrap_position<phase_space>)
                        .def("velocity", &wrap_velocity<phase_space>)
#endif
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_samples_phase_space(lua_State* L)
{
    phase_space<3, float>::luaopen(L);
    phase_space<2, float>::luaopen(L);
    observables::samples::blocking_scheme<phase_space<3, float> >::luaopen(L);
    observables::samples::blocking_scheme<phase_space<2, float> >::luaopen(L);
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
