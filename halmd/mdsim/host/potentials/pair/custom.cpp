/*
 * Copyright © 2008-2013 Felix Höfling
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
#include <halmd/mdsim/host/potentials/pair/custom.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * Initialise custom potential parameters
 */
template <typename float_type>
custom<float_type>::custom(
    matrix_type const& sigma
  , matrix_type const& param2   // FIXME rename param[2-3]
  , matrix_type const& param3
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : sigma_(sigma)
  , param2_(check_shape(param2, sigma))     // FIXME rename param[2-3]
  , param3_(check_shape(param3, sigma))
  , logger_(logger)
{
    // FIXME adjust log messages
    LOG("interaction range: σ = " << sigma_);
    LOG("second potential parameter: p2 = " << param2_);
    LOG("third potential parameter: p3 = " << param3_);
}

template <typename float_type>
void custom<float_type>::luaopen(lua_State* L)
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
                        class_<custom, std::shared_ptr<custom> >("custom")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            // FIXME rename param[2-3]
                            .property("sigma", &custom::sigma)
                            .property("param2", &custom::param2)
                            .property("param3", &custom::param3)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_pair_custom(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    custom<double>::luaopen(L);
    forces::pair_full<3, double, custom<double> >::luaopen(L);
    forces::pair_full<2, double, custom<double> >::luaopen(L);
    truncations::truncations_luaopen<double, custom<double> >(L);
#else
    custom<float>::luaopen(L);
    forces::pair_full<3, float, custom<float> >::luaopen(L);
    forces::pair_full<2, float, custom<float> >::luaopen(L);
    truncations::truncations_luaopen<float, custom<float> >(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class custom<double>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(custom<double>)
#else
template class custom<float>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(custom<float>)
#endif

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::pair::custom<double> >;
template class pair_full<2, double, potentials::pair::custom<double> >;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(double, potentials::pair::custom<double>)
#else
template class pair_full<3, float, potentials::pair::custom<float> >;
template class pair_full<2, float, potentials::pair::custom<float> >;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::custom<float>)
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
