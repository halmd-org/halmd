/*
 * Copyright © 2010 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#include <boost/numeric/ublas/io.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <cmath>
#include <stdexcept>
#include <string>

#include <halmd/mdsim/gpu/potentials/lennard_jones.hpp>
#include <halmd/mdsim/gpu/potentials/lennard_jones_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {

template <typename matrix_type>
static matrix_type const&
check_shape(matrix_type const& m, unsigned int size1, unsigned int size2)
{
    if (m.size1() != size1 || m.size2() != size2) {
        throw std::invalid_argument("parameter matrix has invalid shape");
    }
    return m;
}

/**
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
lennard_jones<float_type>::lennard_jones(
    unsigned int ntype1
  , unsigned int ntype2
  , matrix_type const& cutoff
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , boost::shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(check_shape(epsilon, ntype1, ntype2))
  , sigma_(check_shape(sigma, ntype1, ntype2))
  , r_cut_sigma_(check_shape(cutoff, ntype1, ntype2))
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , sigma2_(element_prod(sigma_, sigma_))
  , en_cut_(ntype1, ntype2)
  , g_param_(ntype1 * ntype2)
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned int i = 0; i < ntype1; ++i) {
        for (unsigned int j = 0; j < ntype2; ++j) {
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
                    class_<lennard_jones, boost::shared_ptr<lennard_jones> >("lennard_jones")
                        .def(constructor<
                            unsigned int
                          , unsigned int
                          , matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , boost::shared_ptr<logger_type>
                        >())
                        .property("r_cut", (matrix_type const& (lennard_jones::*)() const) &lennard_jones::r_cut)
                        .property("r_cut_sigma", &lennard_jones::r_cut_sigma)
                        .property("epsilon", &lennard_jones::epsilon)
                        .property("sigma", &lennard_jones::sigma)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_lennard_jones(lua_State* L)
{
    lennard_jones<float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class lennard_jones<float>;

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd
