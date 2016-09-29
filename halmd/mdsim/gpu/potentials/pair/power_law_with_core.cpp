/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011      Michael Kopp
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
#include <halmd/mdsim/gpu/potentials/pair/local_r4.hpp>
#include <halmd/mdsim/gpu/potentials/pair/power_law_with_core.hpp>
#include <halmd/mdsim/gpu/potentials/pair/power_law_with_core_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

template <typename T, typename S>
static T const&
check_shape(T const& m1, S const& m2)
{
    if (m1.size1() != m2.size1() || m1.size2() != m2.size2()) {
        throw std::invalid_argument("parameter matrix has invalid shape");
    }
    return m1;
}

/**
 * Initialise potential parameters for power law with core
 *
 * @param cutoff cutoff length in units of sigma
 * @param core core radius in units of sigma (potential diverges at core radius)
 * @param epsilon interaction strength in MD units
 * @param sigma interaction range in MD units
 */
template <typename float_type>
power_law_with_core<float_type>::power_law_with_core(
    matrix_type const& cutoff
  , matrix_type const& core
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , uint_matrix_type const& index
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(check_shape(epsilon, epsilon))
  , sigma_(check_shape(sigma, epsilon))
  , index_(check_shape(index, epsilon))
  , r_cut_sigma_(check_shape(cutoff, epsilon))
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , r_core_sigma_(check_shape(core, epsilon))
  , sigma2_(element_prod(sigma_, sigma_))
  , en_cut_(size1(), size2())
  , g_param_(size1() * size2())
  , g_rr_en_cut_(size1() * size2())
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned i = 0; i < en_cut_.size1(); ++i) {
        for (unsigned j = 0; j < en_cut_.size2(); ++j) {
            float_type ri_cut = 1 / r_cut_sigma_(i, j);
            en_cut_(i, j) = epsilon_(i, j) * std::pow(ri_cut, index_(i, j));
        }
    }

    LOG("interaction strength: ε = " << epsilon_);
    LOG("interaction range: σ = " << sigma_);
    LOG("core radius r_core/σ = " << r_core_sigma_);
    LOG("power law index: n = " << index_);
    LOG("cutoff length: r_c/σ = " << r_cut_sigma_);
    LOG("cutoff energy: U = " << en_cut_);

    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p(0); // initialise unused elements as well
        p[power_law_with_core_kernel::EPSILON] = epsilon_.data()[i];
        p[power_law_with_core_kernel::SIGMA2] = sigma2_.data()[i];
        p[power_law_with_core_kernel::CORE_SIGMA] = r_core_sigma_.data()[i];
        p[power_law_with_core_kernel::INDEX] = index_.data()[i];
        param[i] = p;
    }
    cuda::copy(param, g_param_);

    cuda::host::vector<float2> rr_en_cut(g_rr_en_cut_.size());
    for (size_t i = 0; i < rr_en_cut.size(); ++i) {
        rr_en_cut[i].x = rr_cut_.data()[i];
        rr_en_cut[i].y = en_cut_.data()[i];
    }
    cuda::copy(rr_en_cut, g_rr_en_cut_);
}

template <typename float_type>
void power_law_with_core<float_type>::luaopen(lua_State* L)
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
                        class_<power_law_with_core, std::shared_ptr<power_law_with_core> >("power_law_with_core")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , matrix_type const&
                              , matrix_type const&
                              , uint_matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("r_cut", (matrix_type const& (power_law_with_core::*)() const) &power_law_with_core::r_cut)
                            .property("r_cut_sigma", &power_law_with_core::r_cut_sigma)
                            .property("r_core_sigma", &power_law_with_core::r_core_sigma)
                            .property("epsilon", &power_law_with_core::epsilon)
                            .property("sigma", &power_law_with_core::sigma)
                            .property("index", &power_law_with_core::index)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_pair_power_law_with_core(lua_State* L)
{
    power_law_with_core<float>::luaopen(L);
    local_r4<power_law_with_core<float>>::luaopen(L);
    forces::pair_full<3, float, power_law_with_core<float> >::luaopen(L);
    forces::pair_full<2, float, power_law_with_core<float> >::luaopen(L);
    forces::pair_trunc<3, float, power_law_with_core<float> >::luaopen(L);
    forces::pair_trunc<2, float, power_law_with_core<float> >::luaopen(L);
    forces::pair_trunc<3, float, local_r4<power_law_with_core<float> > >::luaopen(L);
    forces::pair_trunc<2, float, local_r4<power_law_with_core<float> > >::luaopen(L);
    return 0;
}

// explicit instantiation
template class power_law_with_core<float>;
template class local_r4<power_law_with_core<float>>;

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
template class pair_full<3, float, potentials::pair::power_law_with_core<float> >;
template class pair_full<2, float, potentials::pair::power_law_with_core<float> >;
template class pair_trunc<3, float, potentials::pair::power_law_with_core<float> >;
template class pair_trunc<2, float, potentials::pair::power_law_with_core<float> >;
template class pair_trunc<3, float, potentials::pair::local_r4<potentials::pair::power_law_with_core<float> > >;
template class pair_trunc<2, float, potentials::pair::local_r4<potentials::pair::power_law_with_core<float> > >;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
