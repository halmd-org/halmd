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
#include <halmd/options.hpp>

namespace halmd { namespace mdsim { namespace host
{

template <unsigned int dimension, typename float_type>
class particle : public mdsim::particle<dimension, float_type>
{
public:
    typedef mdsim::particle<dimension, float_type> _Base;
    typedef vector<float_type, dimension> vector_type;
    typedef std::vector<size_t> neighbor_list;

public:
    particle(options const& vm);
    virtual ~particle() {}

public:
    /** positions, reduced to extended domain box */
    std::vector<vector_type > r;
    /** minimum image vectors */
    std::vector<vector<int, dimension> > image;
    /** velocities */
    std::vector<vector_type > v;
    /** forces */
    std::vector<vector_type > f;
    /** globally unique particle numbers */
    std::vector<unsigned int> tag;
    /** types */
    std::vector<unsigned int> type;
    /** neighbor lists */
    std::vector<neighbor_list> neighbor;

    /** number of particles in simulation box */
    using _Base::nbox;
    /** number of particle types */
    using _Base::ntype;
};

}}} // namespace halmd::mdsim::host

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_HPP */
