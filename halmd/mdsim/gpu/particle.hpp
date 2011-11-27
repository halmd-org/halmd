/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <lua.hpp>
#include <vector>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <unsigned int dimension, typename float_type>
class particle
  : public mdsim::particle<dimension>
{
public:
    typedef mdsim::particle<dimension> _Base;
    typedef typename type_traits<dimension, float_type>::vector_type vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type gpu_vector_type;
    struct defaults;

    static void luaopen(lua_State* L);

    particle(
        std::vector<unsigned int> const& particles
      , std::vector<double> const& mass
      , unsigned int threads = defaults::threads()
    );
    virtual void set();
    void rearrange(cuda::vector<unsigned int> const& g_index);

    /** grid and block dimensions for CUDA calls */
    cuda::config const dim;

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
    /** reverse particle tags */
    cuda::vector<unsigned int> g_reverse_tag;
    /** mass per type */
    cuda::vector<float_type> g_mass;

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

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type rearrange;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

template <unsigned int dimension, typename float_type>
struct particle<dimension, float_type>::defaults
{
    static unsigned int threads();
};

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_HPP */
