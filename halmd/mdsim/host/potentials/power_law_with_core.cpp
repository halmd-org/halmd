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
#include <cmath>
#include <string>

#include <halmd/mdsim/host/potentials/power_law_with_core.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::numeric::ublas;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {

/**
 * Initialise potential parameters
 */
template <typename float_type>
power_law_with_core<float_type>::power_law_with_core(
    unsigned int ntype
  , array<float, 3> const& cutoff
  , array<float, 3> const& core
  , array<float, 3> const& epsilon
  , array<float, 3> const& sigma
  , array<unsigned int, 3> const& index
  , boost::shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(scalar_matrix<float_type>(ntype, ntype, 1))
  , sigma_(scalar_matrix<float_type>(ntype, ntype, 1))
  , index_(ntype, ntype)
  , sigma2_(ntype, ntype)
  , r_cut_(ntype, ntype)
  , r_cut_sigma_(ntype, ntype)
  , rr_cut_(ntype, ntype)
  , r_core_sigma_(ntype, ntype)
  , en_cut_(scalar_matrix<float_type>(ntype, ntype, 0))
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
            sigma2_(i, j) = std::pow(sigma_(i, j), 2);
            r_cut_(i, j) = r_cut_sigma_(i, j) * sigma_(i, j);
            rr_cut_(i, j) = std::pow(r_cut_(i, j), 2);
            // energy shift due to truncation at cutoff length
            en_cut_(i, j) = get<1>((*this)(rr_cut_(i, j), i, j));
        }
    }

    LOG("interaction strength ε = " << epsilon_);
    LOG("interaction range σ = " << sigma_);
    LOG("core radius r_core/σ = " << r_core_sigma_);
    LOG("power law index: n = " << index_);
    LOG("cutoff length: r_c/σ = " << r_cut_sigma_);
    LOG("cutoff energy U = " << en_cut_);
}

template <typename float_type>
void power_law_with_core<float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("potentials")
                [
                    class_<power_law_with_core, boost::shared_ptr<power_law_with_core> >(module_name())
                        .def(constructor<
                            unsigned int                    // ntype
                          , array<float, 3> const&          // cutoff
                          , array<float, 3> const&          // core
                          , array<float, 3> const&          // epsilon
                          , array<float, 3> const&          // sigma
                          , array<unsigned int, 3> const&   // index
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

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_power_law_with_core(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    power_law_with_core<double>::luaopen(L);
#else
    power_law_with_core<float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class power_law_with_core<double>;
#else
template class power_law_with_core<float>;
#endif

} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd
