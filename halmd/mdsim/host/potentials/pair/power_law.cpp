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
#include <stdexcept>
#include <string>

#include <halmd/mdsim/host/forces/pair_full.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/potentials/pair/power_law_hard_core.hpp>
#include <halmd/mdsim/host/potentials/pair/power_law.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/truncations.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * Initialise potential parameters
 */
template <typename float_type>
power_law<float_type>::power_law(
    matrix_type const& epsilon
  , matrix_type const& sigma
  , uint_matrix_type const& index
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , index_(check_shape(index, epsilon))
  , sigma2_(element_prod(sigma_, sigma_))
  , logger_(logger)
{
    LOG("interaction strength ε = " << epsilon_);
    LOG("interaction range σ = " << sigma_);
    LOG("power law index: n = " << index_);
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
                    namespace_("pair")
                    [
                        class_<power_law, std::shared_ptr<power_law> >("power_law")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , uint_matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("epsilon", &power_law::epsilon)
                            .property("sigma", &power_law::sigma)
                            .property("index", &power_law::index)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_pair_power_law(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    power_law<double>::luaopen(L);
    forces::pair_full<3, double, power_law<double> >::luaopen(L);
    forces::pair_full<2, double, power_law<double> >::luaopen(L);
    truncations::truncations_luaopen<double, power_law<double> >(L);


    adapters::hard_core<power_law<double> >::luaopen(L);
    forces::pair_full<3, double, adapters::hard_core<power_law<double> > >::luaopen(L);
    forces::pair_full<2, double, adapters::hard_core<power_law<double> > >::luaopen(L);
    truncations::truncations_luaopen<double, adapters::hard_core<power_law<double> > >(L);
#else
    power_law<float>::luaopen(L);
    forces::pair_full<3, float, power_law<float> >::luaopen(L);
    forces::pair_full<2, float, power_law<float> >::luaopen(L);
    truncations::truncations_luaopen<float, power_law<float> >(L);

    adapters::hard_core<power_law<float> >::luaopen(L);
    forces::pair_full<3, float, adapters::hard_core<power_law<float> > >::luaopen(L);
    forces::pair_full<2, float, adapters::hard_core<power_law<float> > >::luaopen(L);
    truncations::truncations_luaopen<float, adapters::hard_core<power_law<float> > >(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class power_law<double>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(power_law<double>)

template class adapters::hard_core<power_law<double>>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(adapters::hard_core<power_law<double>>)
#else
template class power_law<float>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(power_law<float>)

template class adapters::hard_core<power_law<float>>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(adapters::hard_core<power_law<float>>)
#endif

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::pair::power_law<double>>;
template class pair_full<2, double, potentials::pair::power_law<double>>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(double, potentials::pair::power_law<double>)

template class pair_full<3, double, potentials::pair::adapters::hard_core<potentials::pair::power_law<double>>>;
template class pair_full<2, double, potentials::pair::adapters::hard_core<potentials::pair::power_law<double>>>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(
    double
  , potentials::pair::adapters::hard_core<potentials::pair::power_law<double>>
)
#else
template class pair_full<3, float, potentials::pair::power_law<float>>;
template class pair_full<2, float, potentials::pair::power_law<float>>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float, potentials::pair::power_law<float>)

template class pair_full<3, float, potentials::pair::adapters::hard_core<potentials::pair::power_law<float>>>;
template class pair_full<2, float, potentials::pair::adapters::hard_core<potentials::pair::power_law<float>>>;
HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(
    float
  , potentials::pair::adapters::hard_core<potentials::pair::power_law<float>>
)
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
