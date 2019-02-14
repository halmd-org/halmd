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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_TRUNCATIONS_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_TRUNCATIONS_HPP

#include <boost/preprocessor/seq/for_each.hpp>

#include <halmd/mdsim/gpu/potentials/pair/truncations/force_shifted.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/sharp.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/shifted.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/smooth_r4.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace truncations {

#ifdef USE_GPU_SINGLE_PRECISION
# define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN_SINGLE(r, data, truncation) \
    forces::pair_trunc<3, float, truncation<potential_type> >::luaopen(L);\
    forces::pair_trunc<2, float, truncation<potential_type> >::luaopen(L);
#else
# define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN_SINGLE(r, data, truncation)
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
# define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN_DOUBLE_SINGLE(r, data, truncation) \
    forces::pair_trunc<3, dsfloat, truncation<potential_type> >::luaopen(L);\
    forces::pair_trunc<2, dsfloat, truncation<potential_type> >::luaopen(L);
#else
# define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN_DOUBLE_SINGLE(r, data, truncation)
#endif

#define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN(r, data, truncation) \
    truncation<potential_type>::luaopen(L);\
    _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN_SINGLE(r, data, truncation)\
    _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN_DOUBLE_SINGLE(r, data, truncation)

#define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(r, potential_type, truncation) \
    template class truncations::truncation<potential_type>;

#define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(r, params, truncation) \
    template class pair_trunc<3, BOOST_PP_TUPLE_ELEM(2,0,params), potentials::pair::truncations::truncation<BOOST_PP_TUPLE_ELEM(2,1,params)> >; \
    template class pair_trunc<2, BOOST_PP_TUPLE_ELEM(2,0,params), potentials::pair::truncations::truncation<BOOST_PP_TUPLE_ELEM(2,1,params)> >; \


template<typename potential_type>
void truncations_luaopen(lua_State* L)
{
    BOOST_PP_SEQ_FOR_EACH(_HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_LUAOPEN, _, HALMD_PAIR_POTENTIAL_TRUNCATIONS)
}

#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE(potential_type) \
    BOOST_PP_SEQ_FOR_EACH(_HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE \
                        , potential_type, HALMD_PAIR_POTENTIAL_TRUNCATIONS)

#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES(float_type, potential_type) \
    BOOST_PP_SEQ_FOR_EACH(_HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCES \
                        , (float_type, potential_type), HALMD_PAIR_POTENTIAL_TRUNCATIONS)

} // namespace truncations
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_TRUNCATIONS_HPP */
