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

#ifndef HALMD_MDSIM_HOST_PARTICLE_HPP
#define HALMD_MDSIM_HOST_PARTICLE_HPP

#include <vector>

#include <halmd/math/vector2d.hpp>
#include <halmd/math/vector3d.hpp>
#include <halmd/mdsim/particle.hpp>

namespace halmd { namespace mdsim { namespace host
{

template <unsigned int dimension, typename float_type>
struct particle : mdsim::particle<dimension, float_type>
{
    /** positions, reduced to extended domain box */
    std::vector<vector<float_type, dimension> > r;
    /** minimum image vectors */
    std::vector<vector<int, dimension> > image;
    /** velocities */
    std::vector<vector<float_type, dimension> > v;
    /** forces */
    std::vector<vector<float_type, dimension> > f;
    /** globally unique particle numbers */
    std::vector<unsigned int> tag;
    /** types */
    std::vector<unsigned int> type;
    /** neighbours lists */
    std::vector<std::vector<unsigned int> > neighbour;
};

}}} // namespace halmd::mdsim::host

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_HPP */
