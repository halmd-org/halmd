/*
 * Copyright © 2016 Daniel Kirchner
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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_ADAPTERS_TRUNCATIONS_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_ADAPTERS_TRUNCATIONS_HPP

#include <halmd/mdsim/host/potentials/pair/adapters/force_shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/adapters/sharp.hpp>
#include <halmd/mdsim/host/potentials/pair/adapters/shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/adapters/smooth_r4.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {
namespace adapters {

template<typename float_type, typename potential_type>
void truncations_luaopen(lua_State* L)
{
    smooth_r4<potential_type>::luaopen(L);
    sharp<potential_type>::luaopen(L);
    shifted<potential_type>::luaopen(L);
    force_shifted<potential_type>::luaopen(L);

    forces::pair_trunc<3, float_type, smooth_r4<potential_type> >::luaopen(L);
    forces::pair_trunc<2, float_type, smooth_r4<potential_type> >::luaopen(L);
    forces::pair_trunc<3, float_type, sharp<potential_type> >::luaopen(L);
    forces::pair_trunc<2, float_type, sharp<potential_type> >::luaopen(L);
    forces::pair_trunc<3, float_type, shifted<potential_type> >::luaopen(L);
    forces::pair_trunc<2, float_type, shifted<potential_type> >::luaopen(L);
    forces::pair_trunc<3, float_type, force_shifted<potential_type> >::luaopen(L);
    forces::pair_trunc<2, float_type, force_shifted<potential_type> >::luaopen(L);
}

#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(potential_type) \
    template class adapters::smooth_r4<potential_type>; \
    template class adapters::sharp<potential_type>; \
    template class adapters::shifted<potential_type>; \
    template class adapters::force_shifted<potential_type>

#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float_type, potential_type) \
    template class pair_trunc<3, float_type, potentials::pair::adapters::smooth_r4<potential_type> >; \
    template class pair_trunc<2, float_type, potentials::pair::adapters::smooth_r4<potential_type> >; \
    template class pair_trunc<3, float_type, potentials::pair::adapters::sharp<potential_type> >; \
    template class pair_trunc<2, float_type, potentials::pair::adapters::sharp<potential_type> >; \
    template class pair_trunc<3, float_type, potentials::pair::adapters::shifted<potential_type> >; \
    template class pair_trunc<2, float_type, potentials::pair::adapters::shifted<potential_type> >; \
    template class pair_trunc<3, float_type, potentials::pair::adapters::force_shifted<potential_type> >; \
    template class pair_trunc<2, float_type, potentials::pair::adapters::force_shifted<potential_type> >

} // namespace adapters
} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_ADAPTERS_TRUNCATIONS_HPP */
