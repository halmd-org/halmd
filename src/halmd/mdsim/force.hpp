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

#ifndef HALMD_MDSIM_FORCE_HPP
#define HALMD_MDSIM_FORCE_HPP

#include <vector>

#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/math/vector4d.hpp>

namespace halmd { namespace mdsim
{

template <int dimension, typename float_type>
class force
{
public:
    typedef vector<float_type, dimension> vector_type;
    typedef vector<double, 1 + (dimension - 1) * dimension / 2> virial_type;

public:
    virtual ~force() {};
    virtual void compute() = 0;
    float_type en_pot() { return en_pot_; }
    std::vector<virial_type> const& virial() { return virial_; }

protected:
    /** average potential energy per particle */
    float_type en_pot_;
    /** average virial per particle for each particle type */
    std::vector<virial_type> virial_;
};

}} // namespace halmd::mdsim

#endif /* ! HALMD_MDSIM_FORCE_HPP */
