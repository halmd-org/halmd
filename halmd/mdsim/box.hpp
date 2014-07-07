/*
 * Copyright © 2010-2011 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_BOX_HPP
#define HALMD_MDSIM_BOX_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <lua.hpp>

#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace mdsim {

template <int dimension>
class box
{
private:
    typedef boost::numeric::ublas::zero_matrix<double> zero_matrix_type;

public:
    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef fixed_vector<double, dimension> vector_type;

    /**
     * Construct simulation domain with given edge vectors.
     *
     * http://nongnu.org/h5md/h5md.html#simulation-box
     */
    box(matrix_type const& edges);

    /**
     * Enforce periodic boundary conditions on argument.
     *
     * Assumes that particle position wraps at most once per call.
     *
     * Map coordinates to (-length_half_[i], length_half_[i])
     * which is appropriate too for relative vectors.
     *
     * Return reduction vector in units of box edge lengths.
     *
     * A GPU version is found in halmd/mdsim/gpu/box_kernel.cuh
     */
    template <typename T>
    T reduce_periodic(T& r) const;

    /**
     * Extend periodically reduced distance vector by image vector.
     *
     * This is the inverse of reduce_periodic.
     *
     * A GPU version is found in halmd/mdsim/gpu/box_kernel.cuh
     */
    template <typename T>
    void extend_periodic(T& r, T const& image) const;

    /**
     * Returns edge vectors.
     */
    matrix_type const& edges() const
    {
        return edges_;
    }

    /**
     * Returns coordinates of lowest corner of simulation domain.
     */
    vector_type origin() const;

    /**
     * Returns edge lengths.
     */
    vector_type const& length() const
    {
        return length_;
    }

    /*
     * Calculates volume of box.
     */
    double volume() const;

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** edge vectors of cuboid */
    matrix_type edges_;
    /** edge lengths of cuboid */
    vector_type length_;
    /** store half value for efficient use in reduce_periodic() */
    vector_type length_half_;
};

template <int dimension> template <typename T>
inline T box<dimension>::reduce_periodic(T& r) const
{
    T image;
    for (size_t j = 0; j < dimension; ++j) {
        if (r[j] > length_half_[j]) {
            r[j] -= length_[j];
            image[j] = 1;
        }
        else if (r[j] < -length_half_[j]) {
            r[j] += length_[j];
            image[j] = -1;
        }
        else
            image[j] = 0;
    }
    return image;
}

template <int dimension> template <typename T>
inline void box<dimension>::extend_periodic(T& r, T const& image) const
{
    r += element_prod(image, static_cast<T>(length_));
}

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_BOX_HPP */
