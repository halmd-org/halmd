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
#include <halmd/mdsim/host/potentials/pair/force_shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/morse.hpp>
#include <halmd/mdsim/host/potentials/pair/shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/smooth_r4.hpp>
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
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , r_min_sigma_(check_shape(r_min, epsilon))
  , logger_(logger)
{
    LOG("depth of potential well: ε = " << epsilon_);
    LOG("width of potential well: σ = " << sigma_);
    LOG("position of potential well: r_min / σ = " << r_min_sigma_);
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
                              , std::shared_ptr<logger>
                            >())
                            .property("epsilon", &morse::epsilon)
                            .property("sigma", &morse::sigma)
                            .property("r_min_sigma", &morse::r_min_sigma)
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
    smooth_r4<morse<double>>::luaopen(L);
    shifted<morse<double>>::luaopen(L);
    force_shifted<morse<double>>::luaopen(L);
    forces::pair_full<3, double, morse<double> >::luaopen(L);
    forces::pair_full<2, double, morse<double> >::luaopen(L);
    forces::pair_trunc<3, double, smooth_r4<morse<double> > >::luaopen(L);
    forces::pair_trunc<2, double, smooth_r4<morse<double> > >::luaopen(L);
    forces::pair_trunc<3, double, shifted<morse<double> > >::luaopen(L);
    forces::pair_trunc<2, double, shifted<morse<double> > >::luaopen(L);
    forces::pair_trunc<3, double, force_shifted<morse<double> > >::luaopen(L);
    forces::pair_trunc<2, double, force_shifted<morse<double> > >::luaopen(L);
#else
    morse<float>::luaopen(L);
    smooth_r4<morse<float>>::luaopen(L);
    shifted<morse<float>>::luaopen(L);
    force_shifted<morse<float>>::luaopen(L);
    forces::pair_full<3, float, morse<float> >::luaopen(L);
    forces::pair_full<2, float, morse<float> >::luaopen(L);
    forces::pair_trunc<3, float, smooth_r4<morse<float> > >::luaopen(L);
    forces::pair_trunc<2, float, smooth_r4<morse<float> > >::luaopen(L);
    forces::pair_trunc<3, float, shifted<morse<float> > >::luaopen(L);
    forces::pair_trunc<2, float, shifted<morse<float> > >::luaopen(L);
    forces::pair_trunc<3, float, force_shifted<morse<float> > >::luaopen(L);
    forces::pair_trunc<2, float, force_shifted<morse<float> > >::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class morse<double>;
template class smooth_r4<morse<double>>;
template class shifted<morse<double>>;
template class force_shifted<morse<double>>;
#else
template class morse<float>;
template class smooth_r4<morse<float>>;
template class shifted<morse<float>>;
template class force_shifted<morse<float>>;
#endif

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::pair::morse<double> >;
template class pair_full<2, double, potentials::pair::morse<double> >;
template class pair_trunc<3, double, potentials::pair::smooth_r4<potentials::pair::morse<double> > >;
template class pair_trunc<2, double, potentials::pair::smooth_r4<potentials::pair::morse<double> > >;
template class pair_trunc<3, double, potentials::pair::shifted<potentials::pair::morse<double> > >;
template class pair_trunc<2, double, potentials::pair::shifted<potentials::pair::morse<double> > >;
template class pair_trunc<3, double, potentials::pair::force_shifted<potentials::pair::morse<double> > >;
template class pair_trunc<2, double, potentials::pair::force_shifted<potentials::pair::morse<double> > >;
#else
template class pair_full<3, float, potentials::pair::morse<float> >;
template class pair_full<2, float, potentials::pair::morse<float> >;
template class pair_trunc<3, float, potentials::pair::smooth_r4<potentials::pair::morse<float> > >;
template class pair_trunc<2, float, potentials::pair::smooth_r4<potentials::pair::morse<float> > >;
template class pair_trunc<3, float, potentials::pair::shifted<potentials::pair::morse<float> > >;
template class pair_trunc<2, float, potentials::pair::shifted<potentials::pair::morse<float> > >;
template class pair_trunc<3, float, potentials::pair::force_shifted<potentials::pair::morse<float> > >;
template class pair_trunc<2, float, potentials::pair::force_shifted<potentials::pair::morse<float> > >;
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
