/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2008-2010 Peter Colberg
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

#include <halmd/mdsim/gpu/forces/pair_full.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
#include <halmd/mdsim/gpu/potentials/lennard_jones_linear.hpp>
#include <halmd/mdsim/gpu/potentials/lennard_jones_linear_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {

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
lennard_jones_linear<float_type>::lennard_jones_linear(
    matrix_type const& cutoff
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , std::shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , r_cut_sigma_(check_shape(cutoff, epsilon))
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , sigma2_(element_prod(sigma_, sigma_))
  , en_cut_(size1(), size2())
  , force_cut_(size1(), size2())
  , g_param_(size1() * size2())
  , g_rr_cut_(size1() * size2())
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned int i = 0; i < en_cut_.size1(); ++i) {
        for (unsigned int j = 0; j < en_cut_.size2(); ++j) {
            float_type rri_cut = std::pow(r_cut_sigma_(i, j), -2);
            float_type r6i_cut = rri_cut * rri_cut * rri_cut;
            en_cut_(i, j) = 4 * epsilon_(i, j) * r6i_cut * (r6i_cut - 1);
            force_cut_(i, j) = 24 * epsilon_(i, j) * r6i_cut * (2 * r6i_cut - 1) / r_cut_(i, j);
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential core width: σ = " << sigma_);
    LOG("potential cutoff length: r_c = " << r_cut_sigma_);
    LOG("potential cutoff energy: U_c = " << en_cut_);
    LOG("potential cutoff force: F_c = " << force_cut_);

    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p;
        p[lennard_jones_linear_kernel::EPSILON] = epsilon_.data()[i];
        p[lennard_jones_linear_kernel::SIGMA2] = sigma2_.data()[i];
        p[lennard_jones_linear_kernel::EN_CUT] = en_cut_.data()[i];
        p[lennard_jones_linear_kernel::FORCE_CUT] = force_cut_.data()[i];
        param[i] = p;
    }

    cuda::copy(param, g_param_);
    cuda::copy(rr_cut_.data(), g_rr_cut_);
}

template <typename float_type>
void lennard_jones_linear<float_type>::luaopen(lua_State* L)
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
                    class_<lennard_jones_linear, std::shared_ptr<lennard_jones_linear> >("lennard_jones_linear")
                        .def(constructor<
                            matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , std::shared_ptr<logger_type>
                        >())
                        .property("r_cut", (matrix_type const& (lennard_jones_linear::*)() const) &lennard_jones_linear::r_cut)
                        .property("r_cut_sigma", &lennard_jones_linear::r_cut_sigma)
                        .property("epsilon", &lennard_jones_linear::epsilon)
                        .property("sigma", &lennard_jones_linear::sigma)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_lennard_jones_linear(lua_State* L)
{
    lennard_jones_linear<float>::luaopen(L);
    forces::pair_full<3, float, lennard_jones_linear<float> >::luaopen(L);
    forces::pair_full<2, float, lennard_jones_linear<float> >::luaopen(L);
    forces::pair_trunc<3, float, lennard_jones_linear<float> >::luaopen(L);
    forces::pair_trunc<2, float, lennard_jones_linear<float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class lennard_jones_linear<float>;

} // namespace potentials

namespace forces {

// explicit instantiation of force modules
template class pair_full<3, float, potentials::lennard_jones_linear<float> >;
template class pair_full<2, float, potentials::lennard_jones_linear<float> >;
template class pair_trunc<3, float, potentials::lennard_jones_linear<float> >;
template class pair_trunc<2, float, potentials::lennard_jones_linear<float> >;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
