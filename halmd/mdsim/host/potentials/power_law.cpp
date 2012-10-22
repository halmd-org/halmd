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
#include <cmath>
#include <stdexcept>
#include <string>

#include <halmd/mdsim/forces/trunc/local_r4.hpp>
#include <halmd/mdsim/host/forces/pair_full.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/potentials/power_law.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
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
 * Initialise potential parameters
 */
template <typename float_type>
power_law<float_type>::power_law(
    unsigned int ntype1
  , unsigned int ntype2
  , matrix_type const& cutoff
  , matrix_type const& epsilon
  , matrix_type const& sigma
  , uint_matrix_type const& index
  , std::shared_ptr<logger_type> logger
)
  // allocate potential parameters
  : epsilon_(check_shape(epsilon, ntype1, ntype2))
  , sigma_(check_shape(sigma, ntype1, ntype2))
  , index_(check_shape(index, ntype1, ntype2))
  , sigma2_(element_prod(sigma_, sigma_))
  , r_cut_sigma_(check_shape(cutoff, ntype1, ntype2))
  , r_cut_(element_prod(sigma_, r_cut_sigma_))
  , rr_cut_(element_prod(r_cut_, r_cut_))
  , en_cut_(ntype1, ntype2)
  , logger_(logger)
{
    // energy shift due to truncation at cutoff length
    for (unsigned i = 0; i < en_cut_.size1(); ++i) {
        for (unsigned j = 0; j < en_cut_.size2(); ++j) {
            boost::tie(boost::tuples::ignore, en_cut_(i, j), boost::tuples::ignore) = (*this)(rr_cut_(i, j), i, j);
        }
    }

    LOG("interaction strength ε = " << epsilon_);
    LOG("interaction range σ = " << sigma_);
    LOG("power law index: n = " << index_);
    LOG("cutoff length r_c = " << r_cut_sigma_);
    LOG("cutoff energy U = " << en_cut_);
}

template <typename float_type>
void power_law<float_type>::luaopen(lua_State* L)
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
                    class_<power_law, std::shared_ptr<power_law> >("power_law")
                        .def(constructor<
                            unsigned int
                          , unsigned int
                          , matrix_type const&
                          , matrix_type const&
                          , matrix_type const&
                          , uint_matrix_type const&
                          , std::shared_ptr<logger_type>
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

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_power_law(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    power_law<double>::luaopen(L);
    forces::pair_full<3, double, power_law<double> >::luaopen(L);
    forces::pair_full<2, double, power_law<double> >::luaopen(L);
    forces::pair_trunc<3, double, power_law<double> >::luaopen(L);
    forces::pair_trunc<2, double, power_law<double> >::luaopen(L);
    forces::pair_trunc<3, double, power_law<double>, mdsim::forces::trunc::local_r4<double> >::luaopen(L);
    forces::pair_trunc<2, double, power_law<double>, mdsim::forces::trunc::local_r4<double> >::luaopen(L);
#else
    power_law<float>::luaopen(L);
    forces::pair_full<3, float, power_law<float> >::luaopen(L);
    forces::pair_full<2, float, power_law<float> >::luaopen(L);
    forces::pair_trunc<3, float, power_law<float> >::luaopen(L);
    forces::pair_trunc<2, float, power_law<float> >::luaopen(L);
    forces::pair_trunc<3, float, power_law<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
    forces::pair_trunc<2, float, power_law<float>, mdsim::forces::trunc::local_r4<float> >::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class power_law<double>;
#else
template class power_law<float>;
#endif

} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::power_law<double> >;
template class pair_full<2, double, potentials::power_law<double> >;
template class pair_trunc<3, double, potentials::power_law<double> >;
template class pair_trunc<2, double, potentials::power_law<double> >;
template class pair_trunc<3, double, potentials::power_law<double>, mdsim::forces::trunc::local_r4<double> >;
template class pair_trunc<2, double, potentials::power_law<double>, mdsim::forces::trunc::local_r4<double> >;
#else
template class pair_full<3, float, potentials::power_law<float> >;
template class pair_full<2, float, potentials::power_law<float> >;
template class pair_trunc<3, float, potentials::power_law<float> >;
template class pair_trunc<2, float, potentials::power_law<float> >;
template class pair_trunc<3, float, potentials::power_law<float>, mdsim::forces::trunc::local_r4<float> >;
template class pair_trunc<2, float, potentials::power_law<float>, mdsim::forces::trunc::local_r4<float> >;
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
