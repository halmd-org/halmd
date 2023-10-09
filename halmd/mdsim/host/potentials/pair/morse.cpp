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
#include <string>

#include <halmd/mdsim/host/forces/pair_full.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/potentials/pair/morse.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * Initialise Morse potential parameters
 */
template <typename float_type>
morse<float_type>::morse(
    matrix_type const& epsilon
  , matrix_type const& sigma
  , matrix_type const& r_min
  , matrix_type const& distortion
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , r_min_(check_shape(r_min, epsilon))
  , r_min_sigma_(element_div(r_min, sigma))
  , distortion_(check_shape(distortion, epsilon))
  , logger_(logger)
{
    LOG("depth of potential well: ε = " << epsilon_);
    LOG("width of potential well: σ = " << sigma_);
    LOG("position of potential well: r_min = " << r_min_);
    LOG("distortion factor: B = " << distortion_);
}

template <typename float_type>
void morse<float_type>::luaopen(lua_State* L)
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
                        class_<morse, std::shared_ptr<morse> >("morse")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , matrix_type const&
                              , matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("epsilon", &morse::epsilon)
                            .property("sigma", &morse::sigma)
                            .property("r_min", &morse::r_min)
                            .property("r_min_sigma", &morse::r_min_sigma)
                            .property("distortion", &morse::distortion)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_pair_morse(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    morse<double>::luaopen(L);
    forces::pair_full<3, double, morse<double> >::luaopen(L);
    forces::pair_full<2, double, morse<double> >::luaopen(L);
    truncations::truncations_luaopen<double, morse<double> >(L);
#else
    morse<float>::luaopen(L);
    forces::pair_full<3, float, morse<float> >::luaopen(L);
    forces::pair_full<2, float, morse<float> >::luaopen(L);
    truncations::truncations_luaopen<float, morse<float> >(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class morse<double>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(morse<double>)
#else
template class morse<float>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(morse<float>)
#endif

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::pair::morse<double> >;
template class pair_full<2, double, potentials::pair::morse<double> >;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(double, potentials::pair::morse<double>)
#else
template class pair_full<3, float, potentials::pair::morse<float> >;
template class pair_full<2, float, potentials::pair::morse<float> >;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::morse<float>)
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
