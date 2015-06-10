/*
 * Copyright © 2014 Nicolas Höft
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

#ifndef HALMD_MDSIM_GEOMETRIES_CUBOID_HPP
#define HALMD_MDSIM_GEOMETRIES_CUBOID_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

#ifndef __CUDACC__
# include <lua.hpp>
#endif

namespace halmd {
namespace mdsim {
namespace geometries {

template <int dimension, typename float_type>
class cuboid
{
public:
    typedef fixed_vector<float_type, dimension> vector_type;

#ifndef __CUDACC__
    cuboid(vector_type origin, vector_type length);

    /**
     * Bind class to Lua
     */
    static void luaopen(lua_State* L);
#endif

    /**
     * returns true if the position is within the geometry
     */
    HALMD_GPU_ENABLED bool operator()(vector_type const& r) const;

private:
    vector_type origin_;
    vector_type edge_length_;
};

template<int dimension, typename float_type>
HALMD_GPU_ENABLED bool cuboid<dimension, float_type>::operator()(vector_type const& r) const
{
    bool inside = true;
    vector_type const dr = r - origin_;

    for (int i = 0; i < dimension; ++i) {
        if (dr[i] < 0 || dr[i] > edge_length_[i]) {
            inside = false;
        }
    }

    return inside;
}

} // namespace geometries
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GEOMETRIES_CUBOID_HPP */
