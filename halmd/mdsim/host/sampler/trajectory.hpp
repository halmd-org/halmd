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

#ifndef HALMD_MDSIM_HOST_SAMPLER_TRAJECTORY_HPP
#define HALMD_MDSIM_HOST_SAMPLER_TRAJECTORY_HPP

#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/samples/host/trajectory.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace sampler
{

template <int dimension, typename float_type>
class trajectory
  : public mdsim::samples::host::trajectory<dimension, float_type>
{
public:
    typedef mdsim::samples::host::trajectory<dimension, float_type> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::core<dimension> core_type;

    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;
    boost::shared_ptr<core_type> core;

    static void luaopen(lua_State* L);

    trajectory(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , boost::shared_ptr<core_type> core
    );
    virtual void acquire();

    typedef typename _Base::sample_vector sample_vector;

    /** periodically extended particle positions */
    using _Base::r;
    /** particle velocities */
    using _Base::v;
    /** simulation time when sample was taken */
    using _Base::time;
};

}}} // namespace mdsim::host::sampler

} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_SAMPLER_TRAJECTORY_HPP */
