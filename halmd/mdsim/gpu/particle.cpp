/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <algorithm>
#include <boost/algorithm/string/predicate.hpp>
#include <exception>
#include <numeric>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_kernel.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace halmd::algorithm::gpu; // radix_sort
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * Allocate microscopic system state.
 *
 * @param particles number of particles per type or species
 */
template <unsigned int dimension, typename float_type>
particle<dimension, float_type>::particle(
    vector<unsigned int> const& particles
  , unsigned int threads
)
  : _Base(particles)
  // default CUDA kernel execution dimensions
  , dim(device::validate(cuda::config((nbox + threads - 1) / threads, threads)))
  // allocate global device memory
  , g_r(nbox)
  , g_image(nbox)
  , g_v(nbox)
  , g_f(nbox)
  , g_reverse_tag(nbox)
  // allocate page-locked host memory
  , h_r(nbox)
  , h_image(nbox)
  , h_v(nbox)
  , h_reverse_tag(nbox)
{
    LOG_DEBUG("number of CUDA execution blocks: " << dim.blocks_per_grid());
    LOG_DEBUG("number of CUDA execution threads per block: " << dim.threads_per_block());

    //
    // As the number of threads may exceed the nmber of particles
    // to account for an integer number of threads per block,
    // we need to allocate excess memory for the GPU vectors.
    //
    // The additional memory is allocated using reserve(), which
    // increases the capacity() without changing the size(). The
    // coordinates of these "virtual" particles will be ignored
    // in cuda::copy or cuda::memset calls.
    //
    try {
#ifdef USE_VERLET_DSFUN
        //
        // Double-single precision requires two single precision
        // "words" per coordinate. We use the first part of a GPU
        // vector for the higher (most significant) words of all
        // particle positions or velocities, and the second part for
        // the lower (least significant) words.
        //
        // The additional memory is allocated using reserve(), which
        // increases the capacity() without changing the size().
        //
        // Take care to pass capacity() as an argument to cuda::copy
        // or cuda::memset calls if needed, as the lower words will
        // be ignored in the operation.
        //
        // Particle images remain in single precision as they
        // contain integer values, and otherwise would not matter
        // for the long-time stability of the integrator.
        //
        LOG("integrate using double-single precision");
        g_r.reserve(2 * dim.threads());
        g_v.reserve(2 * dim.threads());
#else
        LOG_WARNING("integrate using single precision");
        g_r.reserve(dim.threads());
        g_v.reserve(dim.threads());
#endif
        g_image.reserve(dim.threads());
        g_f.reserve(dim.threads());
        g_reverse_tag.reserve(dim.threads());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate particles in global device memory");
        throw;
    }

    // initialise 'ghost' particles to zero
    // this avoids potential nonsense computations resulting in denormalised numbers
    cuda::memset(g_r, 0, g_r.capacity());
    cuda::memset(g_v, 0, g_v.capacity());
    cuda::memset(g_f, 0, g_f.capacity());
    cuda::memset(g_image, 0, g_image.capacity());
    cuda::memset(g_reverse_tag, 0, g_reverse_tag.capacity());

    try {
        cuda::copy(nbox, get_particle_kernel<dimension>().nbox);
        cuda::copy(ntype, get_particle_kernel<dimension>().ntype);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy particle parameters to device symbols");
        throw;
    }
}

/**
 * set particle tags and types
 */
template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::set()
{
    try {
        cuda::configure(dim.grid, dim.block);
        cuda::vector<unsigned int> g_ntypes(ntypes.size());
        cuda::copy(ntypes, g_ntypes);
        get_particle_kernel<dimension>().ntypes.bind(g_ntypes);
        get_particle_kernel<dimension>().tag(g_r, g_v);

        cuda::configure(dim.grid, dim.block);
        get_particle_kernel<dimension>().gen_index(g_reverse_tag);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to set particle tags and types");
        throw;
    }
}

/**
 * rearrange particles by permutation
 */
template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::rearrange(cuda::vector<unsigned int> const& g_index)
{
    scoped_timer_type timer(runtime_.rearrange);
    cuda::vector<float4> g_r_buf(nbox);
    cuda::vector<gpu_vector_type> g_image_buf(nbox);
    cuda::vector<float4> g_v_buf(nbox);
    cuda::vector<unsigned int> g_tag(nbox);

    g_r_buf.reserve(g_r.capacity());
    g_image_buf.reserve(g_image.capacity());
    g_v_buf.reserve(g_v.capacity());
    g_tag.reserve(g_reverse_tag.capacity());

    cuda::configure(dim.grid, dim.block);
    get_particle_kernel<dimension>().r.bind(g_r);
    get_particle_kernel<dimension>().image.bind(g_image);
    get_particle_kernel<dimension>().v.bind(g_v);
    get_particle_kernel<dimension>().rearrange(g_index, g_r_buf, g_image_buf, g_v_buf, g_tag);

    g_r_buf.swap(g_r);
    g_image_buf.swap(g_image);
    g_v_buf.swap(g_v);

    radix_sort<unsigned int> sort(nbox, dim.threads_per_block());
    cuda::configure(dim.grid, dim.block);
    get_particle_kernel<dimension>().gen_index(g_reverse_tag);
    sort(g_tag, g_reverse_tag);
}

template <unsigned int dimension, typename float_type>
unsigned int particle<dimension, float_type>::defaults::threads() {
    return 128;
}

template <int dimension, typename float_type>
static int wrap_dimension(particle<dimension, float_type> const&)
{
    return dimension;
}

template <unsigned int dimension, typename float_type>
void particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("particle_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                class_<particle, shared_ptr<_Base>, _Base>(class_name.c_str())
                    .def(constructor<vector<unsigned int> const&>())
                    .def(constructor<vector<unsigned int> const&, unsigned int>())
                    .property("dimension", &wrap_dimension<dimension, float_type>)
                    .scope[
                        namespace_("defaults")
                        [
                            def("threads", &defaults::threads)
                        ]
                      , class_<runtime>("runtime")
                            .def_readonly("rearrange", &runtime::rearrange)
                    ]
                    .def_readonly("runtime", &particle::runtime_)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_particle(lua_State* L)
{
    particle<3, float>::luaopen(L);
    particle<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class particle<3, float>;
template class particle<2, float>;

} // namespace mdsim
} // namespace gpu
} // namespace halmd
