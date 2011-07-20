/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <boost/multi_array.hpp>
#include <numeric>
#include <functional>
#include <lua.hpp>

#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {

template <int dimension>
class box
{
public:
    static void luaopen(lua_State* L);

    typedef typename type_traits<dimension, double>::vector_type vector_type;

    box(
        size_t nbox
      , vector_type const& length
    );
    box(
        size_t nbox
      , double density
      , vector_type const& ratios
    );

    vector_type const& length() const
    {
        return length_;
    }

    double density() const
    {
        return density_;
    }

    double volume() const
    {
        return std::accumulate(length_.begin(), length_.end(), 1., std::multiplies<double>());
    }

    template <typename T>
    T reduce_periodic(T& r) const;

    template <typename T>
    void extend_periodic(T& r, T const& image) const;

    vector_type origin() const
    {
        return -length_half_;
    }

protected:
    /** edge lengths of cuboid */
    vector_type length_;
    /** number density */
    double density_;
    /** store half value for efficient use in reduce_periodic() */
    vector_type length_half_;
};

/**
 * enforce periodic boundary conditions on argument
 *
 * assumes that particle position wraps at most once per call
 *
 * map coordinates to (-length_half_[i], length_half_[i])
 * which is appropriate too for relative vectors
 *
 * return reduction vector in units of box edge lengths
 *
 * A GPU version is found in halmd/mdsim/gpu/box_kernel.cuh
 */
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

/**
 * extend periodically reduced distance vector by image vector
 *
 * This is the inverse of reduce_periodic.
 *
 * A GPU version is found in halmd/mdsim/gpu/box_kernel.cuh
 */
template <int dimension> template <typename T>
inline void box<dimension>::extend_periodic(T& r, T const& image) const
{
    r += element_prod(image, static_cast<T>(length_));
}

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_BOX_HPP */
