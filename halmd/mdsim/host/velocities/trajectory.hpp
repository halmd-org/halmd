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

#ifndef HALMD_MDSIM_HOST_VELOCITIES_TRAJECTORY_HPP
#define HALMD_MDSIM_HOST_VELOCITIES_TRAJECTORY_HPP

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/observables/host/samples/trajectory.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace velocities
{

template <int dimension, typename float_type>
class trajectory
  : public host::velocity<dimension, float_type>
{
public:
    typedef host::velocity<dimension, float_type> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef observables::host::samples::trajectory<dimension, float_type> sample_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<sample_type> sample;

    trajectory(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<sample_type> sample
    );
    virtual void set();
};

}}} // namespace mdsim::host::velocities

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_VELOCITIES_TRAJECTORY_HPP */
