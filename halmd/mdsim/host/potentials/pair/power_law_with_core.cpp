/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011      Michael Kopp
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
#include <stdexcept>
#include <string>

#include <halmd/mdsim/forces/trunc/local_r4.hpp>
#include <halmd/mdsim/host/forces/pair_full.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/potentials/pair/power_law_with_core.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

template <typename T, typename S>
static T const&
check_shape(T const& m1, S const& m2)
{
    if (m1.size1() != m2.size1() || m1.size2() != m2.size2()) {
        throw std::invalid_argument("parameter matrix has invalid shape");
    }
    return m1;
}

/**
 * Initialise potential parameters
 */
template <typename float_type>
power_law_with_core<float_type>::power_law_with_core(
    matrix_type const& cutoff
  , matrix_type const& core
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , uint_matrix_type const& index
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , index_(check_shape(index, epsilon))
  , sigma2_(element_prod(sigma_, sigma_))
  , r_cut_sigma_(check_shape(cutoff, epsilon))
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , r_core_sigma_(check_shape(core, epsilon))
  , en_cut_(size1(), size2())
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned i = 0; i < en_cut_.size1(); ++i) {
        for (unsigned j = 0; j < en_cut_.size2(); ++j) {
            std::tie(std::ignore, en_cut_(i, j)) = (*this)(rr_cut_(i, j), i, j);
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
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("host")
            [
                namespace_("potentials")
                [
                    namespace_("pair")
                    [
                        class_<power_law_with_core, std::shared_ptr<power_law_with_core> >("power_law_with_core")
                            .def(constructor<
                                matrix_type const&          // cutoff
                              , matrix_type const&          // core
                              , matrix_type const&          // epsilon
                              , matrix_type const&          // sigma
                              , uint_matrix_type const&     // index
                              , std::shared_ptr<logger>
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
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_pair_power_law_with_core(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    power_law_with_core<double>::luaopen(L);
    forces::pair_full<3, double, power_law_with_core<double> >::luaopen(L);
    forces::pair_full<2, double, power_law_with_core<double> >::luaopen(L);
    forces::pair_trunc<3, double, power_law_with_core<double> >::luaopen(L);
    forces::pair_trunc<2, double, power_law_with_core<double> >::luaopen(L);
    forces::pair_trunc<3, double, power_law_with_core<double>, mdsim::forces::trunc::local_r4<double> >::luaopen(L);
    forces::pair_trunc<2, double, power_law_with_core<double>, mdsim::forces::trunc::local_r4<double> >::luaopen(L);
#else
    power_law_with_core<float>::luaopen(L);
    forces::pair_full<3, float, power_law_with_core<float> >::luaopen(L);
    forces::pair_full<2, float, power_law_with_core<float> >::luaopen(L);
    forces::pair_trunc<3, float, power_law_with_core<float> >::luaopen(L);
    forces::pair_trunc<2, float, power_law_with_core<float> >::luaopen(L);
    forces::pair_trunc<3, float, power_law_with_core<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
    forces::pair_trunc<2, float, power_law_with_core<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class power_law_with_core<double>;
#else
template class power_law_with_core<float>;
#endif

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::pair::power_law_with_core<double> >;
template class pair_full<2, double, potentials::pair::power_law_with_core<double> >;
template class pair_trunc<3, double, potentials::pair::power_law_with_core<double> >;
template class pair_trunc<2, double, potentials::pair::power_law_with_core<double> >;
template class pair_trunc<3, double, potentials::pair::power_law_with_core<double>, mdsim::forces::trunc::local_r4<double> >;
template class pair_trunc<2, double, potentials::pair::power_law_with_core<double>, mdsim::forces::trunc::local_r4<double> >;
#else
template class pair_full<3, float, potentials::pair::power_law_with_core<float> >;
template class pair_full<2, float, potentials::pair::power_law_with_core<float> >;
template class pair_trunc<3, float, potentials::pair::power_law_with_core<float> >;
template class pair_trunc<2, float, potentials::pair::power_law_with_core<float> >;
template class pair_trunc<3, float, potentials::pair::power_law_with_core<float>, mdsim::forces::trunc::local_r4<float> >;
template class pair_trunc<2, float, potentials::pair::power_law_with_core<float>, mdsim::forces::trunc::local_r4<float> >;
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
