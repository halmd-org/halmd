/*
 * Copyright © 2011  Michael Kopp and Felix Höfling
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

#include <halmd/mdsim/gpu/potentials/power_law_with_core.hpp>
#include <halmd/mdsim/gpu/potentials/power_law_with_core_kernel.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {

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
    unsigned ntype
  , boost::array<float, 3> const& cutoff
  , boost::array<float, 3> const& core
  , boost::array<float, 3> const& epsilon
  , boost::array<float, 3> const& sigma
  , boost::array<unsigned, 3> const& index
  , boost::shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , index_(ntype, ntype)
  , r_cut_sigma_(ntype, ntype)
  , r_cut_(ntype, ntype)
  , rr_cut_(ntype, ntype)
  , r_core_sigma_(ntype, ntype)
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
            index_(i, j) = index[i + j];
            r_cut_sigma_(i, j) = cutoff[i + j];
            r_core_sigma_(i, j) = core[i + j];
        }
    }

    // precalculate derived parameters
    for (unsigned i = 0; i < ntype; ++i) {
        for (unsigned j = i; j < ntype; ++j) {
            r_cut_(i, j) = r_cut_sigma_(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            sigma2_(i, j) = std::pow(sigma_(i, j), 2);
            // energy shift due to truncation at cutoff length
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
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("potentials")
                [
                    class_<power_law_with_core, boost::shared_ptr<power_law_with_core> >(module_name())
                        .def(constructor<
                            unsigned
                          , boost::array<float, 3> const&
                          , boost::array<float, 3> const&
                          , boost::array<float, 3> const&
                          , boost::array<float, 3> const&
                          , boost::array<unsigned, 3> const&
                          , boost::shared_ptr<logger_type>
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
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_potentials_power_law_with_core(lua_State* L)
{
    power_law_with_core<float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class power_law_with_core<float>;

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd
