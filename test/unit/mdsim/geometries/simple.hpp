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

#ifndef TEST_UNIT_MDSIM_GEOMETRY_SIMPLE_HPP
#define TEST_UNIT_MDSIM_GEOMETRY_SIMPLE_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

/**
 * excluded positions if any coordinate is larger
 * than the origin
 */
template<int dimension, typename float_type>
class simple_geometry
{
public:
    typedef halmd::fixed_vector<float_type, dimension> vector_type;
#ifndef __CUDACC__
    simple_geometry(vector_type origin) : origin_(origin) {}
#endif
    HALMD_GPU_ENABLED bool operator()(vector_type const& r) const
    {
        for (int i = 0; i < dimension; ++i) {
            if (r[i] > origin_[i])
                return false;
        }
        return true;
    }
private:
    vector_type origin_;
};

#endif /* ! TEST_UNIT_MDSIM_GEOMETRY_SIMPLE_HPP */
