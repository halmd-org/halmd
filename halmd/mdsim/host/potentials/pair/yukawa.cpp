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
#include <halmd/mdsim/host/potentials/pair/yukawa.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * Initialise Yukawa potential parameters
 */
template <typename float_type>
yukawa<float_type>::yukawa(
    matrix_type const& amplitude
  , matrix_type const& sigma
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : amplitude_(amplitude)
  , sigma_(check_shape(sigma, amplitude))
  , logger_(logger)
{
    LOG("amplitude: A = " << amplitude_);
    LOG("screening length: κ⁻¹ = σ = " << sigma_);
}

template <typename float_type>
void yukawa<float_type>::luaopen(lua_State* L)
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
                        class_<yukawa, std::shared_ptr<yukawa> >("yukawa")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("amplitude", &yukawa::amplitude)
                            .property("sigma", &yukawa::sigma)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_pair_yukawa(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    yukawa<double>::luaopen(L);
    forces::pair_full<3, double, yukawa<double> >::luaopen(L);
    forces::pair_full<2, double, yukawa<double> >::luaopen(L);
    truncations::truncations_luaopen<double, yukawa<double> >(L);
#else
    yukawa<float>::luaopen(L);
    forces::pair_full<3, float, yukawa<float> >::luaopen(L);
    forces::pair_full<2, float, yukawa<float> >::luaopen(L);
    truncations::truncations_luaopen<float, yukawa<float> >(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class yukawa<double>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(yukawa<double>)
#else
template class yukawa<float>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(yukawa<float>)
#endif

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::pair::yukawa<double> >;
template class pair_full<2, double, potentials::pair::yukawa<double> >;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(double, potentials::pair::yukawa<double>)
#else
template class pair_full<3, float, potentials::pair::yukawa<float> >;
template class pair_full<2, float, potentials::pair::yukawa<float> >;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::yukawa<float>)
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
