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

#ifndef HALMD_MDSIM_GPU_SAMPLER_TRAJECTORY_HPP
#define HALMD_MDSIM_GPU_SAMPLER_TRAJECTORY_HPP

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/samples/gpu/trajectory.hpp>
#include <halmd/mdsim/samples/host/trajectory.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace sampler
{

template <typename trajectory_type>
class trajectory;

/**
 * Sample trajectory from GPU memory to GPU memory
 */
template <int dimension, typename float_type>
class trajectory<mdsim::samples::gpu::trajectory<dimension, float_type> >
  : public mdsim::samples::gpu::trajectory<dimension, float_type>
{
public:
    // module definitions
    typedef trajectory _Self;
    typedef mdsim::samples::gpu::trajectory<dimension, float_type> _Base;
    static void depends();
    static void options(po::options_description& desc) {}
    static void select(po::variables_map const& vm) {}

    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::core<dimension> core_type;

    shared_ptr<particle_type> particle;
    shared_ptr<box_type> box;
    shared_ptr<core_type> core;

    trajectory(modules::factory& factory, po::variables_map const& vm);
    virtual ~trajectory() {}
    void acquire();

    typedef typename _Base::sample_vector sample_vector;

    /** periodically extended particle positions */
    using _Base::r;
    /** particle velocities */
    using _Base::v;
    /** simulation time when sample was taken */
    using _Base::time;
};

/**
 * Sample trajectory from GPU memory to host memory
 */
template <int dimension, typename float_type>
class trajectory<mdsim::samples::host::trajectory<dimension, float_type> >
  : public mdsim::samples::host::trajectory<dimension, float_type>
{
public:
    // module definitions
    typedef trajectory _Self;
    typedef mdsim::samples::host::trajectory<dimension, float_type> _Base;
    static void options(po::options_description& desc) {}
    static void depends();
    static void select(po::variables_map const& vm) {}

    typedef gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::core<dimension> core_type;

    shared_ptr<particle_type> particle;
    shared_ptr<box_type> box;
    shared_ptr<core_type> core;

    trajectory(modules::factory& factory, po::variables_map const& vm);
    virtual ~trajectory() {}
    void acquire();

    typedef typename _Base::sample_vector sample_vector;

    /** periodically extended particle positions */
    using _Base::r;
    /** particle velocities */
    using _Base::v;
    /** simulation time when sample was taken */
    using _Base::time;
};

}}} // namespace mdsim::gpu::sampler

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_SAMPLER_TRAJECTORY_HPP */
