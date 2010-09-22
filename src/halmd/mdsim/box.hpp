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

#include <vector>

#include <halmd/utility/module.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim
{

template <int dimension>
class box
{
public:
    // module definitions
    typedef box _Self;
    static void options(po::options_description& desc);
    static void depends();
    static void select(po::options const& vm) {}

    typedef typename type_traits<dimension, double>::vector_type vector_type;
    typedef mdsim::particle<dimension> particle_type;

    shared_ptr<particle_type> particle;

    box(modules::factory& factory, po::options const& vm);
    virtual ~box() {}
    void length(vector_type const& value_type);
    vector_type const& length() { return length_; }
    void density(double value_type);
    double density() { return density_; }

    template <typename T>
    T reduce_periodic(T& r) const;

    vector_type origin() const { return -length_half_; }

protected:
    /** edge lengths of cuboid */
    vector_type length_;
    /** edge lengths of cuboid relative to maximum edge length */
    vector_type scale_;
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
 */
template <int dimension> template <typename T>
inline T box<dimension>::reduce_periodic(T& r) const
{
    vector_type image;
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

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_BOX_HPP */
