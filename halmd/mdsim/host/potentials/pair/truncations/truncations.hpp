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

#ifndef HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_TRUNCATIONS_HPP
#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_TRUNCATIONS_HPP

#include <boost/preprocessor/seq/for_each.hpp>

#include <halmd/mdsim/host/potentials/pair/truncations/force_shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/sharp.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/shifted.hpp>
#include <halmd/mdsim/host/potentials/pair/truncations/smooth_r4.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace potentials {
namespace pair {
namespace truncations {

#define HALMD_PAIR_POTENTIAL_TRUNCATIONS (smooth_r4)(sharp)(shifted)(force_shifted)

#define _HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN(r, data, truncation) \
    truncation<potential_type>::luaopen(L);\
    forces::pair_trunc<3, float_type, truncation<potential_type> >::luaopen(L);\
    forces::pair_trunc<2, float_type, truncation<potential_type> >::luaopen(L);

#define _HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(r, potential_type, truncation) \
    template class truncations::truncation<potential_type>;

#define _HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(r, params, truncation) \
    template class pair_trunc<3, BOOST_PP_TUPLE_ELEM(2,0,params), potentials::pair::truncations::truncation<BOOST_PP_TUPLE_ELEM(2,1,params)> >; \
    template class pair_trunc<2, BOOST_PP_TUPLE_ELEM(2,0,params), potentials::pair::truncations::truncation<BOOST_PP_TUPLE_ELEM(2,1,params)> >; \

template<typename float_type, typename potential_type>
void truncations_luaopen(lua_State* L)
{
    BOOST_PP_SEQ_FOR_EACH(_HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN, _, HALMD_PAIR_POTENTIAL_TRUNCATIONS)
}

#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(potential_type) \
    BOOST_PP_SEQ_FOR_EACH(_HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE \
                        , potential_type, HALMD_PAIR_POTENTIAL_TRUNCATIONS)

#define HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float_type, potential_type) \
    BOOST_PP_SEQ_FOR_EACH(_HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES \
                        , (float_type, potential_type), HALMD_PAIR_POTENTIAL_TRUNCATIONS)

} // namespace truncations
} // namespace pair
} // namespace potentials
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POTENTIALS_PAIR_TRUNCATIONS_TRUNCATIONS_HPP */
