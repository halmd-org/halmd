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

#ifndef HALMD_OBSERVABLES_GPU_TRAJECTORY_HPP
#define HALMD_OBSERVABLES_GPU_TRAJECTORY_HPP

#include <lua.hpp>
#include <map>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/observables/gpu/samples/trajectory.hpp>
#include <halmd/observables/host/samples/trajectory.hpp>
#include <halmd/observables/trajectory.hpp>

namespace halmd
{
namespace observables { namespace gpu
{

template <typename sample_type>
class trajectory;

/**
 * Sample trajectory from GPU memory to GPU memory
 */
template <int dimension, typename float_type>
class trajectory<gpu::samples::trajectory<dimension, float_type> >
  : public observables::trajectory<dimension>
{
public:
    typedef observables::trajectory<dimension> _Base;
    typedef gpu::samples::trajectory<dimension, float_type> sample_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;

    boost::shared_ptr<sample_type> sample;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void luaopen(lua_State* L);

    trajectory(
        boost::shared_ptr<sample_type> sample
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
    );
    virtual void acquire(double time);

    // data stream interface
    virtual void register_request(uint64_t step, boost::function<void(uint64_t)> callback);
    virtual void notify(uint64_t step);

protected:
    // list of data requests from sinks
    typedef std::multimap<uint64_t, boost::function<void(uint64_t)> > request_container_type;
    request_container_type request_;
};

/**
 * Sample trajectory from GPU memory to host memory
 */
template <int dimension, typename float_type>
class trajectory<host::samples::trajectory<dimension, float_type> >
  : public observables::trajectory<dimension>
{
public:
    typedef observables::trajectory<dimension> _Base;
    typedef host::samples::trajectory<dimension, float_type> sample_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef fixed_vector<float_type, dimension> vector_type;

    boost::shared_ptr<sample_type> sample;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<box_type> box;

    static void luaopen(lua_State* L);

    trajectory(
        boost::shared_ptr<sample_type> sample
      , boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
    );
    virtual void acquire(double time);

    // data stream interface
    virtual void register_request(uint64_t step, boost::function<void(uint64_t)> callback);
    virtual void notify(uint64_t step);

protected:
    // list of data requests from sinks
    typedef std::multimap<uint64_t, boost::function<void(uint64_t)> > request_container_type;
    request_container_type request_;
};

}} // namespace observables::gpu

} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_TRAJECTORY_HPP */
