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
#include <halmd/mdsim/gpu/potentials/pair/yukawa.hpp>
#include <halmd/mdsim/gpu/potentials/pair/yukawa_kernel.hpp>
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
yukawa<float_type>::yukawa(
    matrix_type const& amplitude
  , matrix_type const& sigma
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : amplitude_(amplitude)
  , sigma_(check_shape(sigma, amplitude))
  , g_param_(size1() * size2())
  , t_param_(g_param_)
  , logger_(logger)
{
    LOG("amplitude: A = " << amplitude_);
    LOG("screening length: κ⁻¹ = σ = " << sigma_);

    // copy parameters to CUDA device
    cuda::memory::host::vector<float2> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 2> p(0);
        p[yukawa_kernel::AMPLITUDE] = amplitude_.data()[i];
        p[yukawa_kernel::SIGMA] = sigma_.data()[i];
        param[i] = p;
    }
    cuda::copy(param.begin(), param.end(), g_param_.begin());
}

template <typename float_type>
void yukawa<float_type>::luaopen(lua_State* L)
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
                        class_<yukawa, std::shared_ptr<yukawa> >("yukawa")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("amplitude", &yukawa::amplitude)
                            .property("sigma", &yukawa::sigma)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_pair_yukawa(lua_State* L)
{
    yukawa<float>::luaopen(L);
#ifdef USE_GPU_SINGLE_PRECISION
    forces::pair_full<3, float, yukawa<float> >::luaopen(L);
    forces::pair_full<2, float, yukawa<float> >::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    forces::pair_full<3, dsfloat, yukawa<float> >::luaopen(L);
    forces::pair_full<2, dsfloat, yukawa<float> >::luaopen(L);
#endif
    truncations::truncations_luaopen<yukawa<float> >(L);
    return 0;
}

// explicit instantiation
template class yukawa<float>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(yukawa<float>)

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifdef USE_GPU_SINGLE_PRECISION
template class pair_full<3, float, potentials::pair::yukawa<float> >;
template class pair_full<2, float, potentials::pair::yukawa<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::yukawa<float>)
#endif

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class pair_full<3, dsfloat, potentials::pair::yukawa<float> >;
template class pair_full<2, dsfloat, potentials::pair::yukawa<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(dsfloat, potentials::pair::yukawa<float>)
#endif

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
