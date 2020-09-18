/*
 * Copyright © 2008-2013 Felix Höfling
 * Copyright © 2008-2010 Peter Colberg
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
#include <halmd/mdsim/gpu/potentials/pair/morse.hpp>
#include <halmd/mdsim/gpu/potentials/pair/morse_kernel.hpp>
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
morse<float_type>::morse(
    matrix_type const& epsilon
  , matrix_type const& sigma
  , matrix_type const& r_min
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , r_min_sigma_(check_shape(r_min, epsilon))
  , g_param_(size1() * size2())
  , t_param_(g_param_)
  , logger_(logger)
{
    LOG("depth of potential well: ε = " << epsilon_);
    LOG("width of potential well: σ = " << sigma_);
    LOG("position of potential well: r_min / σ = " << r_min_sigma_);

    // copy parameters to CUDA device
    cuda::memory::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p(0);
        p[morse_kernel::EPSILON] = epsilon_.data()[i];
        p[morse_kernel::SIGMA] = sigma_.data()[i];
        p[morse_kernel::R_MIN_SIGMA] = r_min_sigma_.data()[i];
        param[i] = p;
    }
    cuda::copy(param.begin(), param.end(), g_param_.begin());
}

template <typename float_type>
void morse<float_type>::luaopen(lua_State* L)
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
                        class_<morse, std::shared_ptr<morse> >("morse")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("epsilon", &morse::epsilon)
                            .property("sigma", &morse::sigma)
                            .property("r_min_sigma", &morse::r_min_sigma)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_pair_morse(lua_State* L)
{
    morse<float>::luaopen(L);
#ifdef USE_GPU_SINGLE_PRECISION
    forces::pair_full<3, float, morse<float> >::luaopen(L);
    forces::pair_full<2, float, morse<float> >::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    forces::pair_full<3, dsfloat, morse<float> >::luaopen(L);
    forces::pair_full<2, dsfloat, morse<float> >::luaopen(L);
#endif
    truncations::truncations_luaopen<morse<float> >(L);
    return 0;
}

// explicit instantiation
template class morse<float>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(morse<float>)

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifdef USE_GPU_SINGLE_PRECISION
template class pair_full<3, float, potentials::pair::morse<float> >;
template class pair_full<2, float, potentials::pair::morse<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::morse<float>)
#endif

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class pair_full<3, dsfloat, potentials::pair::morse<float> >;
template class pair_full<2, dsfloat, potentials::pair::morse<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(dsfloat, potentials::pair::morse<float>)
#endif

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
