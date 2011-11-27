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

#include <halmd/mdsim/gpu/potentials/morse.hpp>
#include <halmd/mdsim/gpu/potentials/morse_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {

/**
 * Initialise parameters of the potential
 */
template <typename float_type>
morse<float_type>::morse(
    unsigned int ntype1
  , unsigned int ntype2
  , matrix_type const& cutoff
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , matrix_type const& r_min
  , shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(sigma)
  , r_min_sigma_(r_min)
  , r_cut_sigma_(cutoff)
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , en_cut_(ntype1, ntype2)
  , g_param_(ntype1 * ntype2)
  , g_rr_cut_(ntype1 * ntype2)
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned i = 0; i < ntype1; ++i) {
        for (unsigned j = 0; j < ntype2; ++j) {
            float_type a = exp(r_min_sigma_(i, j) - r_cut_sigma_(i, j));
            en_cut_(i, j) = epsilon_(i, j) * (a - 2) * a;
        }
    }

    LOG("depth of potential well: ε = " << epsilon_);
    LOG("width of potential well: σ = " << sigma_);
    LOG("position of potential well: r_min / σ = " << r_min_sigma_);
    LOG("cutoff radius of potential: r_c / σ = " << r_cut_sigma_);
    LOG("potential energy at cutoff: U = " << en_cut_);

    // copy parameters to CUDA device
    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p;
        p[morse_kernel::EPSILON] = epsilon_.data()[i];
        p[morse_kernel::SIGMA] = sigma_.data()[i];
        p[morse_kernel::R_MIN_SIGMA] = r_min_sigma_.data()[i];
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
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("potentials")
                [
                    class_<morse, shared_ptr<morse> >(module_name())
                        .def(constructor<
                            unsigned int
                          , unsigned int
                          , matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , shared_ptr<logger_type>
                        >())
                        .property("r_cut", (matrix_type const& (morse::*)() const) &morse::r_cut)
                        .property("r_cut_sigma", &morse::r_cut_sigma)
                        .property("epsilon", &morse::epsilon)
                        .property("sigma", &morse::sigma)
                        .property("r_min_sigma", &morse::r_min_sigma)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_morse(lua_State* L)
{
    morse<float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class morse<float>;

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd
