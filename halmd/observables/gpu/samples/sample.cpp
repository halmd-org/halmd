/*
 * Copyright 2016 Daniel Kirchner
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HALMD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HALMD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <luaponte/luaponte.hpp>
#include <memory>

#include <halmd/observables/gpu/samples/sample.hpp>
#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace samples {

template <typename sample_type>
static std::size_t wrap_nparticle(sample_type const& self)
{
    return self.data().size();
}

template <typename sample_type>
static std::size_t wrap_dimension(sample_type const&)
{
    return sample_type::dimension;
}

template <typename sample_type>
static std::function<typename sample_type::array_type const&()>
wrap_get(std::shared_ptr<sample_type> self)
{
    return [self]() -> typename sample_type::array_type const& {
        return self->data();
    };
}

template <int dimension, typename data_type>
void sample<dimension, data_type>::luaopen(lua_State* L)
{
    using namespace luaponte;

    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("gpu")
            [
                namespace_("samples")
                [
                    class_<sample, std::shared_ptr<sample> >()
                        .def(constructor<std::size_t>())
                        .property("nparticle", &wrap_nparticle<sample>)
                        .property("dimension", &wrap_dimension<sample>)
                        .def("get", &wrap_get<sample>)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_samples_sample(lua_State* L)
{
    sample<4, float4>::luaopen(L);
    sample<3, float4>::luaopen(L);
    sample<2, float4>::luaopen(L);
    sample<2, float2>::luaopen(L);
    sample<1, float>::luaopen(L);

    observables::samples::blocking_scheme<sample<4, float4> >::luaopen(L);
    observables::samples::blocking_scheme<sample<3, float4> >::luaopen(L);
    observables::samples::blocking_scheme<sample<2, float4> >::luaopen(L);
    observables::samples::blocking_scheme<sample<2, float2> >::luaopen(L);

    sample<4, int4>::luaopen(L);
    sample<3, int4>::luaopen(L);
    sample<2, int4>::luaopen(L);
    sample<2, int2>::luaopen(L);
    sample<1, int>::luaopen(L);

    sample<4, uint4>::luaopen(L);
    sample<3, uint4>::luaopen(L);
    sample<2, uint4>::luaopen(L);
    sample<2, uint2>::luaopen(L);
    sample<1, unsigned int>::luaopen(L);
    return 0;
}


} // namespace samples
} // namespace gpu
} // namespace observables
} // namespace halmd
