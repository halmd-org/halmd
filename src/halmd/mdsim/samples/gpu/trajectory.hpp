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

#ifndef HALMD_MDSIM_SAMPLES_GPU_TRAJECTORY_HPP
#define HALMD_MDSIM_SAMPLES_GPU_TRAJECTORY_HPP

#include <boost/mpl/if.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include <cuda_wrapper.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim { namespace samples { namespace gpu
{

template <int dimension, typename float_type>
class trajectory
{
public:
    // module definitions
    typedef trajectory _Self;
    static void options(po::options_description& desc) {}
    static void depends();
    static void select(po::variables_map const& vm) {}

    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef typename particle_type::gpu_vector_type gpu_vector_type;

    trajectory(modules::factory& factory, po::variables_map const& vm);
    virtual ~trajectory() {}
    virtual void acquire() = 0;

    shared_ptr<particle_type> particle;

    /** sample vector type for all particles of a species */
    typedef cuda::vector<gpu_vector_type> sample_vector;
    /** sample pointer type for all particle of a species */
    typedef shared_ptr<sample_vector> sample_vector_ptr;
    /** sample pointer type for all species */
    typedef std::vector<sample_vector_ptr> sample_vector_ptr_vector;

    /** periodically extended particle positions */
    sample_vector_ptr_vector r;
    /** particle velocities */
    sample_vector_ptr_vector v;
    /** simulation time when sample was taken */
    double time;
};

}}} // namespace mdsim::samples::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_SAMPLES_GPU_TRAJECTORY_HPP */
