/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <halmd/mdsim/module.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/numeric/host/blas/vector.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim
{

template <int dimension>
class box
{
public:
    typedef numeric::host::blas::vector<double, dimension> vector_type;
    typedef mdsim::particle<dimension> particle_type;

public:
    box(options const& vm);
    virtual ~box() {}

    void length(vector_type const& value_type);
    vector_type const& length() { return length_; }
    void density(double value_type);
    double density() { return density_; }
    vector_type reduce_periodic(vector_type& r) const;

public:
    boost::shared_ptr<particle_type> particle;

protected:
    /** edge lengths of cuboid */
    vector_type length_;
    /** store half value for efficient use in reduce_periodic() */
    vector_type length_half_;
    /** edge lengths of cuboid relative to maximum edge length */
    vector_type scale_;
    /** number density */
    double density_;
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
template <int dimension>
inline typename box<dimension>::vector_type
box<dimension>::reduce_periodic(vector_type& r) const
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

template <int dimension>
class module<box<dimension> >
{
public:
    typedef boost::shared_ptr<box<dimension> > pointer;
    static pointer fetch(options const& vm);

private:
    static pointer singleton_;
};

}} // namespace halmd::mdsim

#endif /* ! HALMD_MDSIM_BOX_HPP */
