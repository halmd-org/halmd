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

template <int dimension>
class velocity
{
public:
    typedef vector<double, dimension> vector_type;

public:
    velocity(options const& vm) {}
    virtual ~velocity() {}
    virtual void set() = 0;
};

template <int dimension>
class factory<velocity<dimension> >
{
public:
    typedef boost::shared_ptr<velocity<dimension> > pointer;
    static pointer fetch(options const& vm);

private:
    static pointer singleton_;
};

}} // namespace halmd::mdsim

#endif /* ! HALMD_MDSIM_VELOCITY_HPP */
