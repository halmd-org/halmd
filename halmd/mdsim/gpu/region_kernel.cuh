/*
 * Copyright © 2015 Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_REGION_KERNEL_CUH
#define HALMD_MDSIM_GPU_REGION_KERNEL_CUH

#include <halmd/numeric/blas/fixed_vector.hpp>

template<typename geometry_type>
struct geometry_predicate
{
    typedef typename geometry_type::vector_type vector_type;

    geometry_predicate(float4* position, geometry_type& geometry)
      : position_(position)
      , geometry_(geometry)
    {};

    HALMD_GPU_ENABLED bool operator()(int i)
    {
        vector_type r;
        unsigned int type;
        tie(r, type) <<= position_[i];
        return geometry_(r);
    }

private:
    float4* position_; // position array
    geometry_type const geometry_;

};

#endif /* ! HALMD_MDSIM_GPU_REGION_KERNEL_CUH */
