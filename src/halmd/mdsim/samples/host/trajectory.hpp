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

#ifndef HALMD_MDSIM_SAMPLES_HOST_TRAJECTORY_HPP
#define HALMD_MDSIM_SAMPLES_HOST_TRAJECTORY_HPP

#include <boost/shared_ptr.hpp>
#include <vector>

#include <halmd/mdsim/particle.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>

namespace halmd
{
namespace mdsim { namespace samples { namespace host
{

template <int dimension, typename float_type>
class trajectory
{
public:
    // module definitions
    typedef trajectory _Self;
    static void options(po::options_description& desc) {}
    static void depends();
    static void select(po::options const& vm) {}

    // this module is used by mdsim::gpu::sampler::trajectory
    // and must not depend on host::particle
    typedef mdsim::particle<dimension> particle_type;
    typedef typename mdsim::host::particle<dimension, float_type>::vector_type vector_type;

    trajectory(modules::factory& factory, po::options const& vm);
    virtual ~trajectory() {}
    virtual void acquire() = 0;

    shared_ptr<particle_type> particle;

    /** sample vector type for all particles of a species */
    typedef std::vector<vector_type> sample_vector;
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

}}} // namespace mdsim::samples::host

} // namespace halmd

#endif /* ! HALMD_MDSIM_SAMPLES_HOST_TRAJECTORY_HPP */
