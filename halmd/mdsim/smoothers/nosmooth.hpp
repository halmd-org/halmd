/*
 * Copyright © 2012  Nicolas Höft
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_MDSIM_SMOOTHERS_NOSMOOTH_HPP
#define HALMD_MDSIM_SMOOTHERS_NOSMOOTH_HPP

#include <halmd/config.hpp>

#ifndef __CUDACC__
# include <lua.hpp>
#endif

namespace halmd {
namespace mdsim {
namespace smoothers {

/**
 * This is the smoothing function which does no smoothing at all.
 */
class nosmooth
{
public:
#ifndef __CUDACC__
    static void luaopen(lua_State* L);
#endif

    /**
     * Calculate the smoothing based on the particles distance r,
     * the cutoff distance r_cut, the scalar force |F(r)| f_abs, and
     * the potential U(r).
     * In this case, all variables remain unchanged
     */
    template <typename float_type>
    HALMD_GPU_ENABLED void operator()(float_type, float_type, float_type&, float_type&) const {}
};

} // namespace smoothers
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_SMOOTHERS_NOSMOOTH_HPP */
