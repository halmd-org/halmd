/*
 * Copyright © 2010-2013 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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
#include <halmd/mdsim/host/potentials/pair/force_shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/lennard_jones.hpp>
#include <halmd/mdsim/host/potentials/pair/sharp.hpp>
#include <halmd/mdsim/host/potentials/pair/shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/smooth_r4.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {

/**
 * Initialise Lennard-Jones potential parameters
 */
template <typename float_type>
lennard_jones<float_type>::lennard_jones(
    matrix_type const& epsilon
  , matrix_type const& sigma
  , std::shared_ptr<logger> logger
)
  // allocate potential parameters
  : epsilon_(epsilon)
  , sigma_(check_shape(sigma, epsilon))
  , sigma2_(element_prod(sigma_, sigma_))
  , logger_(logger)
{
    LOG("potential well depths: ε = " << epsilon_);
    LOG("potential core width: σ = " << sigma_);
}

template <typename float_type>
void lennard_jones<float_type>::luaopen(lua_State* L)
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
                        class_<lennard_jones, std::shared_ptr<lennard_jones> >("lennard_jones")
                            .def(constructor<
                                matrix_type const&
                              , matrix_type const&
                              , std::shared_ptr<logger>
                            >())
                            .property("epsilon", &lennard_jones::epsilon)
                            .property("sigma", &lennard_jones::sigma)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_potentials_pair_lennard_jones(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    lennard_jones<double>::luaopen(L);
    smooth_r4<lennard_jones<double>>::luaopen(L);
    sharp<lennard_jones<double>>::luaopen(L);
    shifted<lennard_jones<double>>::luaopen(L);
    force_shifted<lennard_jones<double>>::luaopen(L);
    forces::pair_full<3, double, lennard_jones<double> >::luaopen(L);
    forces::pair_full<2, double, lennard_jones<double> >::luaopen(L);
    forces::pair_trunc<3, double, smooth_r4<lennard_jones<double> > >::luaopen(L);
    forces::pair_trunc<2, double, smooth_r4<lennard_jones<double> > >::luaopen(L);
    forces::pair_trunc<3, double, sharp<lennard_jones<double> > >::luaopen(L);
    forces::pair_trunc<2, double, sharp<lennard_jones<double> > >::luaopen(L);
    forces::pair_trunc<3, double, shifted<lennard_jones<double> > >::luaopen(L);
    forces::pair_trunc<2, double, shifted<lennard_jones<double> > >::luaopen(L);
    forces::pair_trunc<3, double, force_shifted<lennard_jones<double> > >::luaopen(L);
    forces::pair_trunc<2, double, force_shifted<lennard_jones<double> > >::luaopen(L);
#else
    lennard_jones<float>::luaopen(L);
    smooth_r4<lennard_jones<float>>::luaopen(L);
    sharp<lennard_jones<float>>::luaopen(L);
    shifted<lennard_jones<float>>::luaopen(L);
    force_shifted<lennard_jones<float>>::luaopen(L);
    forces::pair_full<3, float, lennard_jones<float> >::luaopen(L);
    forces::pair_full<2, float, lennard_jones<float> >::luaopen(L);
    forces::pair_trunc<3, float, smooth_r4<lennard_jones<float> > >::luaopen(L);
    forces::pair_trunc<2, float, smooth_r4<lennard_jones<float> > >::luaopen(L);
    forces::pair_trunc<3, float, sharp<lennard_jones<float> > >::luaopen(L);
    forces::pair_trunc<2, float, sharp<lennard_jones<float> > >::luaopen(L);
    forces::pair_trunc<3, float, shifted<lennard_jones<float> > >::luaopen(L);
    forces::pair_trunc<2, float, shifted<lennard_jones<float> > >::luaopen(L);
    forces::pair_trunc<3, float, force_shifted<lennard_jones<float> > >::luaopen(L);
    forces::pair_trunc<2, float, force_shifted<lennard_jones<float> > >::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class lennard_jones<double>;
template class smooth_r4<lennard_jones<double>>;
template class sharp<lennard_jones<double>>;
template class shifted<lennard_jones<double>>;
template class force_shifted<lennard_jones<double>>;
#else
template class lennard_jones<float>;
template class smooth_r4<lennard_jones<float>>;
template class sharp<lennard_jones<float>>;
template class shifted<lennard_jones<float>>;
template class force_shifted<lennard_jones<float>>;
#endif

} // namespace pair
} // namespace potentials

namespace forces {

// explicit instantiation of force modules
#ifndef USE_HOST_SINGLE_PRECISION
template class pair_full<3, double, potentials::pair::lennard_jones<double> >;
template class pair_full<2, double, potentials::pair::lennard_jones<double> >;
template class pair_trunc<3, double, potentials::pair::smooth_r4<potentials::pair::lennard_jones<double> > >;
template class pair_trunc<2, double, potentials::pair::smooth_r4<potentials::pair::lennard_jones<double> > >;
template class pair_trunc<3, double, potentials::pair::sharp<potentials::pair::lennard_jones<double> > >;
template class pair_trunc<2, double, potentials::pair::sharp<potentials::pair::lennard_jones<double> > >;
template class pair_trunc<3, double, potentials::pair::shifted<potentials::pair::lennard_jones<double> > >;
template class pair_trunc<2, double, potentials::pair::shifted<potentials::pair::lennard_jones<double> > >;
template class pair_trunc<3, double, potentials::pair::force_shifted<potentials::pair::lennard_jones<double> > >;
template class pair_trunc<2, double, potentials::pair::force_shifted<potentials::pair::lennard_jones<double> > >;
#else
template class pair_full<3, float, potentials::pair::lennard_jones<float> >;
template class pair_full<2, float, potentials::pair::lennard_jones<float> >;
template class pair_trunc<3, float, potentials::pair::smooth_r4<potentials::pair::lennard_jones<float> > >;
template class pair_trunc<2, float, potentials::pair::smooth_r4<potentials::pair::lennard_jones<float> > >;
template class pair_trunc<3, float, potentials::pair::sharp<potentials::pair::lennard_jones<float> > >;
template class pair_trunc<2, float, potentials::pair::sharp<potentials::pair::lennard_jones<float> > >;
template class pair_trunc<3, float, potentials::pair::shifted<potentials::pair::lennard_jones<float> > >;
template class pair_trunc<2, float, potentials::pair::shifted<potentials::pair::lennard_jones<float> > >;
template class pair_trunc<3, float, potentials::pair::force_shifted<potentials::pair::lennard_jones<float> > >;
template class pair_trunc<2, float, potentials::pair::force_shifted<potentials::pair::lennard_jones<float> > >;
#endif

} // namespace forces
} // namespace host
} // namespace mdsim
} // namespace halmd
