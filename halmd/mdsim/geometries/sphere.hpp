/*
 * Copyright © 2019 Roya Ebrahimi Viand
 * Copyright © 2023 Felix Höfling
 * Copyright © 2021 Jaslo Ziska
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
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

#ifndef __CUDACC__
# include <halmd/io/logger.hpp>
# include <lua.hpp>
# include <memory>
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
    sphere(vector_type const& centre, float_type const& radius);

    /**
     * Log geometry information
     *
     * The logger cannot be used as a class member because the host and gpu versions of the class would have different
     * sizes which would lead to memory problems in the gpu code.
     *
     * The logger_ argument must have an underscore in it's name because of the way the LOG() macro works.
     */
    void log(std::shared_ptr<halmd::logger> logger_ = std::make_shared<halmd::logger>()) const;

    /**
     * returns the volume enclosed by the geometry
     */
    float_type volume() const;

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
    vector_type centre_;
    float_type radius_;
    float_type radius2_;
};

template<int dimension, typename float_type>
HALMD_GPU_ENABLED bool sphere<dimension, float_type>::operator()(vector_type const& r) const
{
    bool inside = true;
    vector_type const dr = r - centre_;

    if (inner_prod(dr, dr) > radius2_) {
        inside = false;
    }

    return inside;
}

} // namespace geometries
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GEOMETRIES_SPHERE_HPP */
