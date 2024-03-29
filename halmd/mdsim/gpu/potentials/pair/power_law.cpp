/*
 * Copyright © 2008-2013 Felix Höfling
 * Copyright © 2008-2011 Peter Colberg
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
#include <cmath>
#include <stdexcept>
#include <string>

#include <halmd/mdsim/gpu/forces/pair_full.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/potentials/pair/power_law.hpp>
#include <halmd/mdsim/gpu/potentials/pair/power_law_hard_core.hpp>
#include <halmd/mdsim/gpu/potentials/pair/power_law_kernel.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

/**
 * Initialise power law potential parameters
 */
template <typename float_type>
power_law<float_type>::power_law(
    matrix_type const& epsilon
  , matrix_type const& sigma
  , uint_matrix_type const& index
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , index_(check_shape(index, epsilon))
  , sigma2_(element_prod(sigma_, sigma_))
  , g_param_(size1() * size2())
  , t_param_(g_param_)
  , logger_(logger)
{
    LOG("interaction strength: ε = " << epsilon_);
    LOG("interaction range: σ = " << sigma_);
    LOG("power law index: n = " << index_);

    cuda::memory::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p(0); // initialise unused elements as well
        p[power_law_kernel::EPSILON] = epsilon_.data()[i];
        p[power_law_kernel::SIGMA2] = sigma2_.data()[i];
        p[power_law_kernel::INDEX] = index_.data()[i];
        param[i] = p;
    }
    cuda::copy(param.begin(), param.end(), g_param_.begin());
}

template <typename float_type>
void power_law<float_type>::luaopen(lua_State* L)
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
                        class_<power_law, std::shared_ptr<power_law> >("power_law")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , uint_matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("epsilon", &power_law::epsilon)
                            .property("sigma", &power_law::sigma)
                            .property("index", &power_law::index)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_pair_power_law(lua_State* L)
{
    power_law<float>::luaopen(L);
#ifdef USE_GPU_SINGLE_PRECISION
    forces::pair_full<3, float, power_law<float> >::luaopen(L);
    forces::pair_full<2, float, power_law<float> >::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    forces::pair_full<3, dsfloat, power_law<float> >::luaopen(L);
    forces::pair_full<2, dsfloat, power_law<float> >::luaopen(L);
#endif
    truncations::truncations_luaopen<power_law<float> >(L);

    adapters::hard_core<power_law<float>>::luaopen(L);
#ifdef USE_GPU_SINGLE_PRECISION
    forces::pair_full<3, float, adapters::hard_core<power_law<float> > >::luaopen(L);
    forces::pair_full<2, float, adapters::hard_core<power_law<float> > >::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    forces::pair_full<3, dsfloat, adapters::hard_core<power_law<float> > >::luaopen(L);
    forces::pair_full<2, dsfloat, adapters::hard_core<power_law<float> > >::luaopen(L);
#endif
    truncations::truncations_luaopen<adapters::hard_core<power_law<float> > >(L);

    return 0;
}

// explicit instantiation
template class power_law<float>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(power_law<float>)

template class adapters::hard_core<power_law<float>>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(adapters::hard_core<power_law<float>>)

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifdef USE_GPU_SINGLE_PRECISION
template class pair_full<3, float, potentials::pair::power_law<float> >;
template class pair_full<2, float, potentials::pair::power_law<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::power_law<float>)

template class pair_full<3, float, potentials::pair::adapters::hard_core<potentials::pair::power_law<float> > >;
template class pair_full<2, float, potentials::pair::adapters::hard_core<potentials::pair::power_law<float> > >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(
    float
  , potentials::pair::adapters::hard_core<potentials::pair::power_law<float> >
  )
#endif

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class pair_full<3, dsfloat, potentials::pair::power_law<float> >;
template class pair_full<2, dsfloat, potentials::pair::power_law<float> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(dsfloat, potentials::pair::power_law<float>)

template class pair_full<3, dsfloat, potentials::pair::adapters::hard_core<potentials::pair::power_law<float> > >;
template class pair_full<2, dsfloat, potentials::pair::adapters::hard_core<potentials::pair::power_law<float> > >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(
  dsfloat
  , potentials::pair::adapters::hard_core<potentials::pair::power_law<float> >
)
#endif

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
