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

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <cmath>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/forces/morse.hpp>
#include <halmd/mdsim/gpu/forces/morse_kernel.hpp>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{

/**
 * Initialise parameters of the potential
 */
template <typename float_type>
morse<float_type>::morse(
    unsigned ntype
  , array<float, 3> const& cutoff
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
  , array<float, 3> const& r_min
)
  // allocate potential parameters
  : epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , r_min_(ntype, ntype)
  , en_cut_(ntype, ntype)
  , r_cut_(ntype, ntype)
  , rr_cut_(ntype, ntype)
  , g_param_(epsilon_.data().size())
  , g_rr_cut_(epsilon_.data().size())
{
    // FIXME support any number of types
    matrix_type r_cut_sigma(ntype, ntype);
    matrix_type r_min_sigma(ntype, ntype);
    for (unsigned i = 0; i < std::min(ntype, 2U); ++i) {
        for (unsigned j = i; j < std::min(ntype, 2U); ++j) {
            epsilon_(i, j) = epsilon[i + j];
            sigma_(i, j) = sigma[i + j];
            r_cut_sigma(i, j) = cutoff[i + j];
            r_min_sigma(i, j) = r_min[i + j];
        }
    }

    // precalculate derived parameters
    for (unsigned i = 0; i < ntype; ++i) {
        for (unsigned j = i; j < ntype; ++j) {
            r_min_(i, j) = r_min_sigma(i, j) * sigma_(i, j);
            r_cut_(i, j) = r_cut_sigma(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            // energy shift due to truncation at cutoff length
            float_type a = exp((r_min_(i, j) * r_min_(i, j) - rr_cut_(i, j)) / sigma_(i, j));
            en_cut_(i, j) = epsilon_(i, j) * (a - 2) * a;
        }
    }

    LOG("depth of potential well: ε = " << epsilon_);
    LOG("width of potential well: σ = " << sigma_);
    LOG("position of potential well: r_min / σ = " << r_min_sigma);
    LOG("cutoff radius of potential: r_c / σ = " << r_cut_sigma);
    LOG("potential energy at cutoff: U = " << en_cut_);

    // copy parameters to CUDA device
    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p;
        p[morse_kernel::EPSILON] = epsilon_.data()[i];
        p[morse_kernel::SIGMA] = sigma_.data()[i];
        p[morse_kernel::R_MIN] = r_min_.data()[i];
        p[morse_kernel::EN_CUT] = en_cut_.data()[i];
        param[i] = p;
    }
    cuda::copy(param, g_param_);
    cuda::copy(rr_cut_.data(), g_rr_cut_);
}

template <typename float_type>
void morse<float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "halmd_wrapper")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("forces")
                [
                    class_<morse, shared_ptr<morse> >(module_name())
                        .def(constructor<
                            unsigned
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , array<float, 3> const&
                          , array<float, 3> const&
                        >())
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(0) //< distance of derived to base class
    [
        &morse<float>::luaopen
    ];

    lua_wrapper::register_(2) //< distance of derived to base class
    [
        &pair_trunc<3, float, morse<float> >::luaopen
    ]
    [
        &pair_trunc<2, float, morse<float> >::luaopen
    ];
}

// explicit instantiation
template class morse<float>;
template class pair_trunc<3, float, morse<float> >;
template class pair_trunc<2, float, morse<float> >;

}}} // namespace mdsim::gpu::forces

} // namespace halmd
