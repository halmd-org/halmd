/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#include <halmd/mdsim/forces/trunc/local_r4.hpp>
#include <halmd/mdsim/gpu/forces/pair_full.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones.hpp>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {

template <typename matrix_type>
static matrix_type const&
check_shape(matrix_type const& m1, matrix_type const& m2)
{
    if (m1.size1() != m2.size1() || m1.size2() != m2.size2()) {
        throw std::invalid_argument("parameter matrix has invalid shape");
    }
    return m1;
}

/**
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
lennard_jones<float_type>::lennard_jones(
    matrix_type const& cutoff
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , r_cut_sigma_(check_shape(cutoff, epsilon))
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , sigma2_(element_prod(sigma_, sigma_))
  , en_cut_(size1(), size2())
  , g_param_(size1() * size2())
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned int i = 0; i < en_cut_.size1(); ++i) {
        for (unsigned int j = 0; j < en_cut_.size2(); ++j) {
            float_type rri_cut = std::pow(r_cut_sigma_(i, j), -2);
            float_type r6i_cut = rri_cut * rri_cut * rri_cut;
            en_cut_(i, j) = 4 * epsilon_(i, j) * r6i_cut * (r6i_cut - 1);
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential core width: σ = " << sigma_);
    LOG("potential cutoff length: r_c = " << r_cut_sigma_);
    LOG("potential cutoff energy: U = " << en_cut_);

    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p;
        p[lennard_jones_kernel::EPSILON] = epsilon_.data()[i];
        p[lennard_jones_kernel::RR_CUT] = rr_cut_.data()[i];
        p[lennard_jones_kernel::SIGMA2] = sigma2_.data()[i];
        p[lennard_jones_kernel::EN_CUT] = en_cut_.data()[i];
        param[i] = p;
    }

    cuda::copy(param, g_param_);
}

template <typename float_type>
void lennard_jones<float_type>::luaopen(lua_State* L)
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
                        class_<lennard_jones, std::shared_ptr<lennard_jones> >("lennard_jones")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("r_cut", (matrix_type const& (lennard_jones::*)() const) &lennard_jones::r_cut)
                            .property("r_cut_sigma", &lennard_jones::r_cut_sigma)
                            .property("epsilon", &lennard_jones::epsilon)
                            .property("sigma", &lennard_jones::sigma)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_pair_lennard_jones(lua_State* L)
{
    lennard_jones<float>::luaopen(L);
    forces::pair_full<3, float, lennard_jones<float> >::luaopen(L);
    forces::pair_full<2, float, lennard_jones<float> >::luaopen(L);
    forces::pair_trunc<3, float, lennard_jones<float> >::luaopen(L);
    forces::pair_trunc<2, float, lennard_jones<float> >::luaopen(L);
    forces::pair_trunc<3, float, lennard_jones<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
    forces::pair_trunc<2, float, lennard_jones<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class lennard_jones<float>;

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
template class pair_full<3, float, potentials::pair::lennard_jones<float> >;
template class pair_full<2, float, potentials::pair::lennard_jones<float> >;
template class pair_trunc<3, float, potentials::pair::lennard_jones<float> >;
template class pair_trunc<2, float, potentials::pair::lennard_jones<float> >;
template class pair_trunc<3, float, potentials::pair::lennard_jones<float>, mdsim::forces::trunc::local_r4<float> >;
template class pair_trunc<2, float, potentials::pair::lennard_jones<float>, mdsim::forces::trunc::local_r4<float> >;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
