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

#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/factory.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
class box
{
public:
    typedef vector<float_type, dimension> vector_type;
    typedef particle<dimension, float_type> particle_type;

public:
    box(options const& vm);
    virtual ~box() {}

    void length(vector_type const& value_type);
    vector_type const& length() { return length_; }
    void density(float_type value_type);
    float_type density() { return density_; }

public:
    boost::shared_ptr<particle_type> particle;

protected:
    /** edge lengths of cuboid */
    vector_type length_;
    /** edge lengths of cuboid relative to maximum edge length */
    vector_type scale_;
    /** number density */
    float_type density_;
};

template <int dimension, typename float_type>
class factory<box<dimension, float_type> >
{
public:
    typedef boost::shared_ptr<box<dimension, float_type> > box_ptr;
    static box_ptr fetch(options const& vm);

private:
    static box_ptr box_;
};

}} // namespace halmd::mdsim

#endif /* ! HALMD_MDSIM_BOX_HPP */
