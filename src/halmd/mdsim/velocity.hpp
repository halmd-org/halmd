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

#ifndef HALMD_MDSIM_VELOCITY_HPP
#define HALMD_MDSIM_VELOCITY_HPP

#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/factory.hpp>
#include <halmd/options.hpp>

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
class velocity
{
public:
    typedef vector<float_type, dimension> vector_type;

public:
    velocity(options const& vm) {}
    virtual ~velocity() {}
    virtual void set() = 0;
};

template <int dimension, typename float_type>
class factory<velocity<dimension, float_type> >
{
public:
    typedef boost::shared_ptr<velocity<dimension, float_type> > velocity_ptr;
    static velocity_ptr fetch(options const& vm);

private:
    static velocity_ptr velocity_;
};

}} // namespace halmd::mdsim

#endif /* ! HALMD_MDSIM_VELOCITY_HPP */
