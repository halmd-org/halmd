/*
 * Copyright © 2019 Roya Ebrahimi Viand
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

#ifndef HALMD_MDSIM_GEOMETRIES_SPHERE_HPP
#define HALMD_MDSIM_GEOMETRIES_SPHERE_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
// #include <halmd/numeric/blas/blas1.hpp>

#ifndef __CUDACC__
# include <halmd/io/logger.hpp>  // FIXME remove if not needed
# include <lua.hpp>
#endif

namespace halmd {
namespace mdsim {
namespace geometries {

template <int dimension, typename float_type>
class sphere
{
public:
    typedef fixed_vector<float_type, dimension> vector_type;

#ifndef __CUDACC__
   
    sphere(vector_type center, float_type radius);
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
    vector_type center_;
    float_type radius_;
};

template<int dimension, typename float_type>
HALMD_GPU_ENABLED bool sphere<dimension, float_type>::operator()(vector_type const& r) const
{
    bool inside = true;
    vector_type const dr = r - center_;

#ifndef __CUDACC__
    LOG_TRACE("geometry: dr = " << dr);
#endif
    
    float_type sum=0;
    float_type r2= radius_*radius_;

    for (int i = 0; i < dimension; ++i) {
	sum+=dr[i]*dr[i];
        }
    if ( sum > r2) {
	inside = false;
    }
   /* if ( inner_prod(dr, dr) > r2 ) {
	inside = false;
	}*/
   /*for (int i = 0; i < dimension; ++i) {
        if ( abs(dr[i]) > radius_) {
            inside = false;
        }
    }*/

    return inside;
}

} // namespace geometries
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GEOMETRIES_SPHERE_HPP */
