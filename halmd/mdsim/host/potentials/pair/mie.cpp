/*
 * Copyright © 2008-2023 Felix Höfling
 * Copyright © 2008-2011 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/forces/pair_full.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/potentials/pair/mie.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost::numeric::ublas;
using namespace std;

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
mie<float_type>::mie(
    matrix_type const& epsilon
  , matrix_type const& sigma
  , uint_matrix_type const& index_m
  , uint_matrix_type const& index_n
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , epsilon_C_(epsilon)                              // multiply prefactor C below
  , sigma_(check_shape(sigma, epsilon))
  , index_m_(check_shape(index_m, epsilon))
  , index_m_2_(index_m_ / 2)
  , index_n_(check_shape(index_n, epsilon))
  , index_n_2_(index_n_ / 2)
  , sigma2_(element_prod(sigma_, sigma_))
  , logger_(logger)
{
    LOG("potential well depths: ε = " << epsilon_);
    LOG("interaction range: σ = " << sigma_);
    LOG("index of repulsion: m = " << index_m_);
    LOG("index of attraction: n = " << index_n_);


    // check conditions on power law indices (after logging output)
    for (unsigned i = 0; i < index_m_.size1(); ++i) {
        for (unsigned j = 0; j < index_m_.size2(); ++j) {
            // indices must be even
            if (index_m_(i, j) & 1 || index_n_(i, j) & 1) {
                throw std::logic_error("power law indices of potential must be even");
            }
            if (index_m_(i, j) <= index_n_(i, j)) {
                throw std::logic_error("repulsive part of potential must be stronger than attraction");
            }
        }
    }

    // compute prefactor C(m,n) and multiply with epsilon
    for (unsigned i = 0; i < index_m_.size1(); ++i) {
        for (unsigned j = 0; j < index_m_.size2(); ++j) {
            float_type m = index_m_(i, j);  // promote to floating-point numbers
            float_type n = index_n_(i, j);
            epsilon_C_(i, j) *= m / (m - n) * std::pow(m / n, n / (m - n));
        }
    }
}

template <typename float_type>
void mie<float_type>::luaopen(lua_State* L)
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
                        class_<mie, std::shared_ptr<mie> >("mie")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , uint_matrix_type const&
                              , uint_matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("epsilon", &mie::epsilon)
                            .property("sigma", &mie::sigma)
                            .property("index_m", &mie::index_m)
                            .property("index_n", &mie::index_n)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_pair_mie(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    mie<double>::luaopen(L);
    forces::pair_full<3, double, mie<double> >::luaopen(L);
    forces::pair_full<2, double, mie<double> >::luaopen(L);
    truncations::truncations_luaopen<double, mie<double> >(L);
#else
    mie<float>::luaopen(L);
    forces::pair_full<3, float, mie<float> >::luaopen(L);
    forces::pair_full<2, float, mie<float> >::luaopen(L);
    truncations::truncations_luaopen<float, mie<float> >(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class mie<double>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(mie<double>)
#else
template class mie<float>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(mie<float>)
#endif

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::pair::mie<double> >;
template class pair_full<2, double, potentials::pair::mie<double> >;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(double, potentials::pair::mie<double>)
#else
template class pair_full<3, float, potentials::pair::mie<float> >;
template class pair_full<2, float, potentials::pair::mie<float> >;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::mie<float>)
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
