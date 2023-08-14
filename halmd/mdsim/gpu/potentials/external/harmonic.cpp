/*
 * Copyright © 2014 Felix Höfling
 * Copyright © 2020 Jaslo Ziska
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

#include <boost/numeric/ublas/io.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <stdexcept>
#include <string>

#include <halmd/mdsim/gpu/forces/external.hpp>
#include <halmd/mdsim/gpu/potentials/external/harmonic.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {

/**
 * Initialise harmonic potential
 */
template <int dimension, typename float_type>
harmonic<dimension, float_type>::harmonic(
    scalar_container_type const& stiffness
  , vector_container_type const& offset
  , shared_ptr<logger> logger
)
  // allocate potential parameters
  : stiffness_(stiffness)
  , offset_(offset)
  , g_param_(stiffness.size())
  , t_param_(g_param_)
  , logger_(logger)
{
    // check parameter size
    if (stiffness_.size() != offset_.size()) {
        throw invalid_argument("parameter lists have mismatching shapes");
    }

    LOG("potential stiffness: K = " << stiffness_);
    LOG("potential offset: r₀ = " << offset_);

    // merge parameters in a single array and copy to device
    cuda::memory::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        param[i] <<= tie(offset_[i], stiffness_[i]);
    }
    cuda::copy(param.begin(), param.end(), g_param_.begin());
}

template <int dimension, typename float_type>
void harmonic<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static string class_name("harmonic_" + to_string(dimension));
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("potentials")
                [
                    namespace_("external")
                    [
                        class_<harmonic, shared_ptr<harmonic>>(class_name.c_str())
                            .def(constructor<
                                scalar_container_type const&
                              , vector_container_type const&
                              , shared_ptr<logger>
                             >())
                            .property("stiffness", &harmonic::stiffness)
                            .property("offset", &harmonic::offset)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_external_harmonic(lua_State* L)
{
    harmonic<3, float>::luaopen(L);
    harmonic<2, float>::luaopen(L);
#ifdef USE_GPU_SINGLE_PRECISION
    forces::external<3, float, harmonic<3, float>>::luaopen(L);
    forces::external<2, float, harmonic<2, float>>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    forces::external<3, dsfloat, harmonic<3, float>>::luaopen(L);
    forces::external<2, dsfloat, harmonic<2, float>>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
template class harmonic<3, float>;
template class harmonic<2, float>;

} // namespace external
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
using namespace potentials::external;

#ifdef USE_GPU_SINGLE_PRECISION
template class external<3, float, harmonic<3, float>>;
template class external<2, float, harmonic<2, float>>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class external<3, dsfloat, harmonic<3, float>>;
template class external<2, dsfloat, harmonic<2, float>>;
#endif

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
