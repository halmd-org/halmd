/*
 * Copyright © 2008-2012  Felix Höfling
 * Copyright © 2008-2010  Peter Colberg
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
#include <string>

#include <halmd/mdsim/gpu/potentials/lennard_jones_linear.hpp>
#include <halmd/mdsim/gpu/potentials/lennard_jones_linear_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {

/**
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
lennard_jones_linear<float_type>::lennard_jones_linear(
    unsigned ntype
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
  , shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , r_cut_sigma_(ntype, ntype)
  , r_cut_(ntype, ntype)
  , rr_cut_(ntype, ntype)
  , sigma2_(ntype, ntype)
  , en_cut_(ntype, ntype)
  , force_cut_(ntype, ntype)
  , g_param_(epsilon_.data().size())
  , g_rr_cut_(epsilon_.data().size())
  , logger_(logger)
{
    // FIXME support any number of types
    for (unsigned i = 0; i < std::min(ntype, 2U); ++i) {
        for (unsigned j = i; j < std::min(ntype, 2U); ++j) {
            epsilon_(i, j) = epsilon[i + j];
            sigma_(i, j) = sigma[i + j];
            r_cut_sigma_(i, j) = cutoff[i + j];
        }
    }

    // precalculate derived parameters
    for (unsigned i = 0; i < ntype; ++i) {
        for (unsigned j = i; j < ntype; ++j) {
            r_cut_(i, j) = r_cut_sigma_(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            sigma2_(i, j) = std::pow(sigma_(i, j), 2);
            // energy shift due to truncation at cutoff length
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
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("potentials")
                [
                    class_<lennard_jones_linear, shared_ptr<lennard_jones_linear> >(module_name())
                        .def(constructor<
                            unsigned
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , shared_ptr<logger_type>
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
    return 0;
}

// explicit instantiation
template class lennard_jones_linear<float>;

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd
