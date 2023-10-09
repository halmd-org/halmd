/*
 * Copyright © 2023 Felix Höfling
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

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <cmath>
#include <string>

#include <halmd/mdsim/gpu/forces/pair_full.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/potentials/pair/custom.hpp>
#include <halmd/mdsim/gpu/potentials/pair/custom_kernel.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * Initialise parameters of the potential
 */
template <typename float_type>
custom<float_type>::custom(
    matrix_type const& sigma
  , matrix_type const& param2
  , matrix_type const& param3
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : sigma_(sigma)
  , param2_(check_shape(param2, sigma))
  , param3_(check_shape(param3, sigma))
  , g_param_(size1() * size2())
  , t_param_(g_param_)
  , logger_(logger)
{
    // FIXME adjust log messages
    LOG("interaction range: σ = " << sigma_);
    LOG("second potential parameter: p2 = " << param2_);
    LOG("third potential parameter: p3 = " << param3_);

    // copy parameters to CUDA device
    cuda::memory::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p(0);
        p[custom_kernel::SIGMA] = sigma_.data()[i];
        p[custom_kernel::PARAM2] = param2_.data()[i];
        p[custom_kernel::PARAM3] = param3_.data()[i];
        param[i] = p;
    }
    cuda::copy(param.begin(), param.end(), g_param_.begin());
}

template <typename float_type>
void custom<float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("potentials")
                [
                    namespace_("pair")
                    [
                        class_<custom, std::shared_ptr<custom> >("custom")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            // FIXME rename param[1-3]
                            .property("sigma", &custom::sigma)
                            .property("param2", &custom::param2)
                            .property("param3", &custom::param3)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_pair_custom(lua_State* L)
{
    custom<float>::luaopen(L);
#ifdef USE_GPU_SINGLE_PRECISION
    forces::pair_full<3, float, custom<float> >::luaopen(L);
    forces::pair_full<2, float, custom<float> >::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    forces::pair_full<3, dsfloat, custom<float> >::luaopen(L);
    forces::pair_full<2, dsfloat, custom<float> >::luaopen(L);
#endif
    truncations::truncations_luaopen<custom<float> >(L);
    return 0;
}

// explicit instantiation
template class custom<float>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(custom<float>)

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifdef USE_GPU_SINGLE_PRECISION
template class pair_full<3, float, potentials::pair::custom<float> >;
template class pair_full<2, float, potentials::pair::custom<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::custom<float>)
#endif

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class pair_full<3, dsfloat, potentials::pair::custom<float> >;
template class pair_full<2, dsfloat, potentials::pair::custom<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(dsfloat, potentials::pair::custom<float>)
#endif

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
