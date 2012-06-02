/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/potentials/power_law.hpp>
#include <halmd/mdsim/gpu/potentials/power_law_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

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
 * Initialise power law potential parameters
 */
template <typename float_type>
power_law<float_type>::power_law(
    unsigned ntype1
  , unsigned ntype2
  , matrix_type const& cutoff
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , uint_matrix_type const& index
  , boost::shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(check_shape(epsilon, ntype1, ntype2))
  , sigma_(check_shape(sigma, ntype1, ntype2))
  , index_(check_shape(index, ntype1, ntype2))
  , r_cut_sigma_(check_shape(cutoff, ntype1, ntype2))
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , sigma2_(element_prod(sigma_, sigma_))
  , en_cut_(ntype1, ntype2)
  , g_param_(ntype1 * ntype2)
  , g_rr_en_cut_(ntype1 * ntype2)
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned i = 0; i < ntype1; ++i) {
        for (unsigned j = 0; j < ntype2; ++j) {
            float_type ri_cut = 1 / r_cut_sigma_(i, j);
            en_cut_(i, j) = epsilon_(i, j) * std::pow(ri_cut, index_(i, j));
        }
    }

    LOG("interaction strength: ε = " << epsilon_);
    LOG("interaction range: σ = " << sigma_);
    LOG("power law index: n = " << index_);
    LOG("cutoff length: r_c = " << r_cut_sigma_);
    LOG("cutoff energy: U = " << en_cut_);

    cuda::host::vector<float4> param(g_param_.size());
    for (size_t i = 0; i < param.size(); ++i) {
        fixed_vector<float, 4> p(0); // initialise unused elements as well
        p[power_law_kernel::EPSILON] = epsilon_.data()[i];
        p[power_law_kernel::SIGMA2] = sigma2_.data()[i];
        p[power_law_kernel::INDEX] = index_.data()[i];
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
                    class_<power_law, boost::shared_ptr<power_law> >(module_name())
                        .def(constructor<
                            unsigned
                          , unsigned
                          , matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , uint_matrix_type const&
                          , boost::shared_ptr<logger_type>
                        >())
                        .property("r_cut", (matrix_type const& (power_law::*)() const) &power_law::r_cut)
                        .property("r_cut_sigma", &power_law::r_cut_sigma)
                        .property("epsilon", &power_law::epsilon)
                        .property("sigma", &power_law::sigma)
                        .property("index", &power_law::index)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_power_law(lua_State* L)
{
    power_law<float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class power_law<float>;

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd
