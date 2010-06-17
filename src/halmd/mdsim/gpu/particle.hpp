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

#ifndef HALMD_MDSIM_GPU_PARTICLE_HPP
#define HALMD_MDSIM_GPU_PARTICLE_HPP

#include <boost/mpl/if.hpp>
#include <vector>

#include <cuda_wrapper.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/numeric/host/blas/vector.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/utility/gpu/device.hpp>

namespace halmd
{
namespace mdsim { namespace gpu
{

template <unsigned int dimension, typename float_type>
class particle
  : public mdsim::particle<dimension>
{
public:
    // module definitions
    typedef particle _Self;
    typedef mdsim::particle<dimension> _Base;
    static void depends();
    static void select(po::options const& vm);
    static void options(po::options_description& desc) {}

    typedef typename boost::mpl::if_c<dimension == 3, float4, float2>::type gpu_vector_type;
    typedef numeric::host::blas::vector<float_type, dimension> vector_type;
    typedef utility::gpu::device device_type;

    shared_ptr<device_type> device;

    particle(po::options const& vm);
    virtual ~particle() {}

    //
    // particles in global device memory
    //

    /** positions, types */
    cuda::vector<float4> g_r;
    /** minimum image vectors */
    cuda::vector<gpu_vector_type> g_image;
    /** velocities, tags */
    cuda::vector<float4> g_v;
    /** forces */
    cuda::vector<gpu_vector_type> g_f;
    /** neighbour lists */
    cuda::vector<unsigned int> g_neighbour;

    //
    // particles in page-locked host memory
    //

    /** positions, types */
    cuda::host::vector<float4> h_r;
    /** minimum image vectors */
    cuda::host::vector<gpu_vector_type> h_image;
    /** velocities, tags */
    cuda::host::vector<float4> h_v;

    /** number of particles in simulation box */
    using _Base::nbox;
    /** number of particle types */
    using _Base::ntype;
    /** number of particles per type */
    using _Base::ntypes;
    /** number of placeholders per neighbour list */
    unsigned int neighbour_size;
    /** neighbour list stride */
    unsigned int neighbour_stride;
    /** grid and block dimensions for CUDA calls */
    const cuda::config dim;
};

}} // namespace mdsim::gpu

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_HPP */
