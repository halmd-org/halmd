/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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
#include <halmd/utility/signal.hpp>

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
template <int dimension, typename float_type>
particle<dimension, float_type>::particle(size_t nparticle, unsigned int threads)
  // default CUDA kernel execution dimensions
  : dim(device::validate(cuda::config((nparticle + threads - 1) / threads, threads)))
  // allocate global device memory
  , g_position_(nparticle)
  , g_image_(nparticle)
  , g_velocity_(nparticle)
  , g_reverse_tag_(nparticle)
  , g_tag_(nparticle)
  , g_force_(nparticle)
  , g_en_pot_(nparticle)
  , g_stress_pot_(nparticle)
  , g_hypervirial_(nparticle)
  // disable auxiliary variables by default
  , aux_flag_(false)
  , aux_valid_(false)
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
        g_position_.reserve(2 * dim.threads());
        g_velocity_.reserve(2 * dim.threads());
#else
        LOG_WARNING("integrate using single precision");
        g_position_.reserve(dim.threads());
        g_velocity_.reserve(dim.threads());
#endif
        g_image_.reserve(dim.threads());
        g_tag_.reserve(dim.threads());
        g_reverse_tag_.reserve(dim.threads());
        g_force_.reserve(dim.threads());
        g_en_pot_.reserve(dim.threads());
        g_stress_pot_.reserve(dim.threads());
        g_hypervirial_.reserve(dim.threads());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate particles in global device memory");
        throw;
    }

    // initialise 'ghost' particles to zero
    // this avoids potential nonsense computations resulting in denormalised numbers
    cuda::memset(g_position_, 0, g_position_.capacity());
    cuda::memset(g_velocity_, 0, g_velocity_.capacity());
    cuda::memset(g_image_, 0, g_image_.capacity());
    cuda::memset(g_tag_, 0, g_tag_.capacity());
    cuda::memset(g_reverse_tag_, 0, g_reverse_tag_.capacity());

    // set particle masses to unit mass
    set_mass(1);

    try {
        cuda::copy(g_tag_.size(), get_particle_kernel<dimension>().nbox);
        cuda::copy(nspecies_, get_particle_kernel<dimension>().ntype);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy particle parameters to device symbols");
        throw;
    }

    LOG("number of particles: " << g_tag_.size());
    LOG("number of particle placeholders: " << g_tag_.capacity());
    LOG("number of particle species: " << nspecies_);
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::aux_enable()
{
    LOG_TRACE("enable computation of auxiliary variables");
    aux_flag_ = true;
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::set_mass(float_type mass)
{
    cuda::configure(dim.grid, dim.block);
    get_particle_kernel<dimension>().set_mass(g_velocity_, g_velocity_.size(), mass);
    cuda::thread::synchronize();
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::prepare()
{
    LOG_TRACE("zero forces");
    cuda::memset(g_force_, 0, g_force_.capacity());

    // indicate whether auxiliary variables are computed this step
    aux_valid_ = aux_flag_;

    if (aux_flag_) {
        LOG_TRACE("zero auxiliary variables");
        cuda::memset(g_en_pot_, 0, g_en_pot_.capacity());
        cuda::memset(g_stress_pot_, 0, g_stress_pot_.capacity());
        cuda::memset(g_hypervirial_, 0, g_hypervirial_.capacity());
        aux_flag_ = false;
    }
}

/**
 * set particle types
 */
template <int dimension, typename float_type>
void particle<dimension, float_type>::set()
{
    try {
        cuda::configure(dim.grid, dim.block);
        get_particle_kernel<dimension>().tag(g_position_);

        cuda::configure(dim.grid, dim.block);
        get_particle_kernel<dimension>().gen_index(g_reverse_tag_);
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
template <int dimension, typename float_type>
void particle<dimension, float_type>::rearrange(cuda::vector<unsigned int> const& g_index)
{
    scoped_timer_type timer(runtime_.rearrange);
    cuda::vector<float4> g_r_buf(g_tag_.size());
    cuda::vector<gpu_vector_type> g_image_buf(g_tag_.size());
    cuda::vector<float4> g_v_buf(g_tag_.size());
    cuda::vector<unsigned int> g_tag(g_tag_.size());

    g_r_buf.reserve(g_position_.capacity());
    g_image_buf.reserve(g_image_.capacity());
    g_v_buf.reserve(g_velocity_.capacity());
    g_tag.reserve(g_reverse_tag_.capacity());

    cuda::configure(dim.grid, dim.block);
    get_particle_kernel<dimension>().r.bind(g_position_);
    get_particle_kernel<dimension>().image.bind(g_image_);
    get_particle_kernel<dimension>().v.bind(g_velocity_);
    get_particle_kernel<dimension>().rearrange(g_index, g_r_buf, g_image_buf, g_v_buf, g_tag);

    g_r_buf.swap(g_position_);
    g_image_buf.swap(g_image_);
    g_v_buf.swap(g_velocity_);

    radix_sort<unsigned int> sort(g_tag_.size(), dim.threads_per_block());
    cuda::configure(dim.grid, dim.block);
    get_particle_kernel<dimension>().gen_index(g_reverse_tag_);
    sort(g_tag, g_reverse_tag_);
}

template <int dimension, typename float_type>
unsigned int particle<dimension, float_type>::defaults::threads() {
    return 128;
}

template <int dimension, typename float_type>
static int wrap_dimension(particle<dimension, float_type> const&)
{
    return dimension;
}

template <typename particle_type>
static typename signal<void ()>::slot_function_type
wrap_aux_enable(shared_ptr<particle_type> self)
{
    return bind(&particle_type::aux_enable, self);
}

template <typename particle_type>
static typename signal<void ()>::slot_function_type
wrap_prepare(shared_ptr<particle_type> self)
{
    return bind(&particle_type::prepare, self);
}

template <typename particle_type>
typename signal<void ()>::slot_function_type
wrap_set(shared_ptr<particle_type> particle)
{
    return bind(&particle_type::set, particle);
}

template <typename particle_type>
struct wrap_particle
  : particle_type
  , luabind::wrap_base
{
    wrap_particle(size_t nparticle, unsigned int threads) : particle_type(nparticle, threads) {}

    wrap_particle(size_t nparticle) : particle_type(nparticle) {}
};

template <int dimension, typename float_type>
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
                class_<particle, shared_ptr<particle>, wrap_particle<particle> >(class_name.c_str())
                    .def(constructor<size_t, unsigned int>())
                    .def(constructor<size_t>())
                    .property("nparticle", &particle::nparticle)
                    .property("nspecies", &particle::nspecies)
                    .def("set_mass", &particle::set_mass)
                    .property("dimension", &wrap_dimension<dimension, float_type>)
                    .property("aux_enable", &wrap_aux_enable<particle>)
                    .property("prepare", &wrap_prepare<particle>)
                    .property("set", &wrap_set<particle>)
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
