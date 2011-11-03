/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/potentials/modified_lennard_jones.hpp>
#include <halmd/mdsim/gpu/potentials/modified_lennard_jones_kernel.hpp>
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
modified_lennard_jones<float_type>::modified_lennard_jones(
    unsigned ntype
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
  , array<unsigned, 3> const& index_m
  , array<unsigned, 3> const& index_n
  , shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , index_m_(ntype, ntype)
  , index_n_(ntype, ntype)
  , r_cut_sigma_(ntype, ntype)
  , r_cut_(ntype, ntype)
  , rr_cut_(ntype, ntype)
  , sigma2_(ntype, ntype)
  , en_cut_(ntype, ntype)
  , g_param_(epsilon_.data().size())
  , g_rr_en_cut_(epsilon_.data().size())
  , logger_(logger)
{
    // FIXME support any number of types
    for (unsigned i = 0; i < std::min(ntype, 2U); ++i) {
        for (unsigned j = i; j < std::min(ntype, 2U); ++j) {
            epsilon_(i, j) = epsilon[i + j];
            sigma_(i, j) = sigma[i + j];
            index_m_(i, j) = index_m[i + j];
            index_n_(i, j) = index_n[i + j];
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
            unsigned m_2 = index_m_(i, j) / 2;
            unsigned n_2 = index_n_(i, j) / 2;
            float_type rni_cut = std::pow(rri_cut, n_2);
            float_type rmni_cut = std::pow(rri_cut, m_2 - n_2);
            en_cut_(i, j) = 4 * epsilon_(i, j) * rni_cut * (rmni_cut - 1);
        }
    }

    LOG("potential well depths: ε = " << epsilon_);
    LOG("interaction range: σ = " << sigma_);
    LOG("index of repulsion: m = " << index_m_);
    LOG("index of attraction: n = " << index_n_);
    LOG("cutoff length: r_c = " << r_cut_sigma_);
    LOG("cutoff energy: U = " << en_cut_);

    // check conditions on power low indices (after logging output)
    for (unsigned i = 0; i < ntype; ++i) {
        for (unsigned j = i; j < ntype; ++j) {
            // indices must be even
            if (index_m_(i, j) & 1 || index_n_(i, j) & 1) {
                throw std::logic_error("power law indices of potential must be even");
            }
            if (index_m_(i, j) <= index_n_(i, j)) {
                throw std::logic_error("repulsive part of potential must be stronger than attraction");
            }
        }
    }

    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p;
        p[modified_lennard_jones_kernel::EPSILON] = epsilon_.data()[i];
        p[modified_lennard_jones_kernel::SIGMA2] = sigma2_.data()[i];
        p[modified_lennard_jones_kernel::INDEX_M_2] = index_m_.data()[i] / 2;
        p[modified_lennard_jones_kernel::INDEX_N_2] = index_n_.data()[i] / 2;
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
void modified_lennard_jones<float_type>::luaopen(lua_State* L)
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
                    class_<modified_lennard_jones, shared_ptr<modified_lennard_jones> >(module_name())
                        .def(constructor<
                            unsigned
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , array<unsigned, 3> const&
                          , array<unsigned, 3> const&
                          , shared_ptr<logger_type>
                        >())
                        .property("r_cut", (matrix_type const& (modified_lennard_jones::*)() const) &modified_lennard_jones::r_cut)
                        .property("r_cut_sigma", &modified_lennard_jones::r_cut_sigma)
                        .property("epsilon", &modified_lennard_jones::epsilon)
                        .property("sigma", &modified_lennard_jones::sigma)
                        .property("index_m", &modified_lennard_jones::index_m)
                        .property("index_n", &modified_lennard_jones::index_n)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_modified_lennard_jones(lua_State* L)
{
    modified_lennard_jones<float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class modified_lennard_jones<float>;

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd
