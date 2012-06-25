/*
 * Copyright © 2008-2012 Peter Colberg
 * Copyright © 2010-2012 Felix Höfling
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

#include <halmd/config.hpp>

#include <halmd/algorithm/gpu/iota.hpp>
#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_kernel.hpp>
#include <halmd/utility/gpu/device.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/signal.hpp>

#include <boost/function_output_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <luaponte/out_value_policy.hpp>

#include <algorithm>
#include <exception>
#include <iterator>
#include <numeric>

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * Allocate microscopic system state.
 *
 * @param particles number of particles per type or species
 */
template <int dimension, typename float_type>
particle<dimension, float_type>::particle(size_type nparticle, unsigned int nspecies)
  // FIXME default CUDA kernel execution dimensions
  : dim(device::validate(cuda::config((nparticle + 128 - 1) / 128, 128)))
  // allocate global device memory
  , nparticle_(nparticle)
  , nspecies_(std::max(nspecies, 1u))
  , g_position_(nparticle)
  , g_image_(nparticle)
  , g_velocity_(nparticle)
  , g_tag_(nparticle)
  , g_reverse_tag_(nparticle)
{
    cache_proxy<position_array_type> g_position = g_position_;
    cache_proxy<image_array_type> g_image = g_image_;
    cache_proxy<velocity_array_type> g_velocity = g_velocity_;
    cache_proxy<tag_array_type> g_tag = g_tag_;
    cache_proxy<reverse_tag_array_type> g_reverse_tag = g_reverse_tag_;

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
        g_position->reserve(2 * dim.threads());
        g_velocity->reserve(2 * dim.threads());
#else
        LOG_WARNING("integrate using single precision");
        g_position->reserve(dim.threads());
        g_velocity->reserve(dim.threads());
#endif
        g_image->reserve(dim.threads());
        g_tag->reserve(dim.threads());
        g_reverse_tag->reserve(dim.threads());
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to allocate particles in global device memory");
        throw;
    }

    // initialise 'ghost' particles to zero
    // this avoids potential nonsense computations resulting in denormalised numbers
    cuda::memset(g_position->begin(), g_position->begin() + g_position->capacity(), 0);
    cuda::memset(g_velocity->begin(), g_velocity->begin() + g_velocity->capacity(), 0);
    cuda::memset(g_image->begin(), g_image->begin() + g_image->capacity(), 0);
    iota(g_tag->begin(), g_tag->begin() + g_tag->capacity(), 0);
    iota(g_reverse_tag->begin(), g_reverse_tag->begin() + g_reverse_tag->capacity(), 0);

    // set particle masses to unit mass
    set_mass(
        *this
      , boost::make_transform_iterator(boost::counting_iterator<tag_type>(0), [](tag_type) {
            return 1;
        })
    );

    try {
        cuda::copy(nparticle_, get_particle_kernel<dimension>().nbox);
        cuda::copy(nspecies_, get_particle_kernel<dimension>().ntype);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy particle parameters to device symbols");
        throw;
    }

    LOG("number of particles: " << nparticle_);
    LOG("number of particle placeholders: " << dim.threads());
    LOG("number of particle species: " << nspecies_);
}

/**
 * rearrange particles by permutation
 */
template <int dimension, typename float_type>
void particle<dimension, float_type>::rearrange(cuda::vector<unsigned int> const& g_index)
{
    cache_proxy<position_array_type> g_position = g_position_;
    cache_proxy<image_array_type> g_image = g_image_;
    cache_proxy<velocity_array_type> g_velocity = g_velocity_;
    cache_proxy<tag_array_type> g_tag = g_tag_;
    cache_proxy<reverse_tag_array_type> g_reverse_tag = g_reverse_tag_;

    scoped_timer_type timer(runtime_.rearrange);

    cuda::vector<float4> position(nparticle_);
    cuda::vector<gpu_vector_type> image(nparticle_);
    cuda::vector<float4> velocity(nparticle_);
    cuda::vector<unsigned int> tag(nparticle_);

    position.reserve(g_position->capacity());
    image.reserve(g_image->capacity());
    velocity.reserve(g_velocity->capacity());
    tag.reserve(g_reverse_tag->capacity());

    cuda::configure(dim.grid, dim.block);
    get_particle_kernel<dimension>().r.bind(*g_position);
    get_particle_kernel<dimension>().image.bind(*g_image);
    get_particle_kernel<dimension>().v.bind(*g_velocity);
    get_particle_kernel<dimension>().tag.bind(*g_tag);
    get_particle_kernel<dimension>().rearrange(g_index, position, image, velocity, tag);

    position.swap(*g_position);
    image.swap(*g_image);
    velocity.swap(*g_velocity);
    cuda::copy(tag.begin(), tag.begin() + tag.capacity(), g_tag->begin());

    iota(g_reverse_tag->begin(), g_reverse_tag->begin() + g_reverse_tag->capacity(), 0);
    radix_sort(tag.begin(), tag.end(), g_reverse_tag->begin());
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::position_type> ()>
wrap_get_position(std::shared_ptr<particle_type const> self)
{
    return [=]() -> std::vector<typename particle_type::position_type> {
        std::vector<typename particle_type::position_type> output;
        output.reserve(self->nparticle());
        get_position(*self, back_inserter(output));
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<void (std::vector<typename particle_type::position_type> const&)>
wrap_set_position(std::shared_ptr<particle_type> self)
{
    return [=](std::vector<typename particle_type::position_type> const& input) {
        if (input.size() != self->nparticle()) {
            throw std::invalid_argument("input array size not equal to number of particles");
        }
        set_position(*self, input.begin());
    };
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::image_type> ()>
wrap_get_image(std::shared_ptr<particle_type const> self)
{
    return [=]() -> std::vector<typename particle_type::image_type> {
        std::vector<typename particle_type::image_type> output;
        output.reserve(self->nparticle());
        get_image(*self, back_inserter(output));
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<void (std::vector<typename particle_type::image_type> const&)>
wrap_set_image(std::shared_ptr<particle_type> self)
{
    return [=](std::vector<typename particle_type::image_type> const& input) {
        if (input.size() != self->nparticle()) {
            throw std::invalid_argument("input array size not equal to number of particles");
        }
        set_image(*self, input.begin());
    };
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::velocity_type> ()>
wrap_get_velocity(std::shared_ptr<particle_type const> self)
{
    return [=]() -> std::vector<typename particle_type::velocity_type> {
        std::vector<typename particle_type::velocity_type> output;
        output.reserve(self->nparticle());
        get_velocity(*self, back_inserter(output));
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<void (std::vector<typename particle_type::velocity_type> const&)>
wrap_set_velocity(std::shared_ptr<particle_type> self)
{
    return [=](std::vector<typename particle_type::velocity_type> const& input) {
        if (input.size() != self->nparticle()) {
            throw std::invalid_argument("input array size not equal to number of particles");
        }
        set_velocity(*self, input.begin());
    };
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::tag_type> ()>
wrap_get_tag(std::shared_ptr<particle_type const> self)
{
    return [=]() -> std::vector<typename particle_type::tag_type> {
        std::vector<typename particle_type::tag_type> output;
        output.reserve(self->nparticle());
        get_tag(
            *self
          , boost::make_function_output_iterator([&](typename particle_type::tag_type t) {
                output.push_back(t + 1);
            })
        );
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<void (std::vector<typename particle_type::tag_type> const&)>
wrap_set_tag(std::shared_ptr<particle_type> self)
{
    typedef typename particle_type::tag_type tag_type;
    return [=](std::vector<tag_type> const& input) {
        if (input.size() != self->nparticle()) {
            throw std::invalid_argument("input array size not equal to number of particles");
        }
        tag_type nparticle = self->nparticle();
        set_tag(
            *self
          , boost::make_transform_iterator(input.begin(), [&](tag_type t) -> tag_type {
                if (t < 1 || t > nparticle) {
                    throw std::invalid_argument("invalid particle tag");
                }
                return t - 1;
            })
        );
    };
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::reverse_tag_type> ()>
wrap_get_reverse_tag(std::shared_ptr<particle_type const> self)
{
    typedef typename particle_type::reverse_tag_type reverse_tag_type;
    return [=]() -> std::vector<reverse_tag_type> {
        std::vector<typename particle_type::reverse_tag_type> output;
        output.reserve(self->nparticle());
        get_reverse_tag(
            *self
          , boost::make_function_output_iterator([&](reverse_tag_type i) {
                output.push_back(i + 1);
            })
        );
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<void (std::vector<typename particle_type::reverse_tag_type> const&)>
wrap_set_reverse_tag(std::shared_ptr<particle_type> self)
{
    typedef typename particle_type::reverse_tag_type reverse_tag_type;
    return [=](std::vector<reverse_tag_type> const& input) {
        if (input.size() != self->nparticle()) {
            throw std::invalid_argument("input array size not equal to number of particles");
        }
        reverse_tag_type nparticle = self->nparticle();
        set_reverse_tag(
            *self
          , boost::make_transform_iterator(input.begin(), [&](reverse_tag_type i) -> reverse_tag_type {
                if (i < 1 || i > nparticle) {
                    throw std::invalid_argument("invalid particle reverse tag");
                }
                return i - 1;
            })
        );
    };
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::species_type> ()>
wrap_get_species(std::shared_ptr<particle_type const> self)
{
    typedef typename particle_type::species_type species_type;
    return [=]() -> std::vector<species_type> {
        std::vector<species_type> output;
        output.reserve(self->nparticle());
        get_species(
            *self
          , boost::make_function_output_iterator([&](typename particle_type::species_type s) {
                output.push_back(s + 1);
            })
        );
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<void (std::vector<typename particle_type::species_type> const&)>
wrap_set_species(std::shared_ptr<particle_type> self)
{
    typedef typename particle_type::species_type species_type;
    return [=](std::vector<species_type> const& input) {
        if (input.size() != self->nparticle()) {
            throw std::invalid_argument("input array size not equal to number of particles");
        }
        species_type nspecies = self->nspecies();
        set_species(
            *self
          , boost::make_transform_iterator(input.begin(), [&](species_type s) -> species_type {
                if (s < 1 || s > nspecies) {
                    throw std::invalid_argument("invalid particle species");
                }
                return s - 1;
            })
        );
    };
}

template <typename particle_type>
static std::function<std::vector<typename particle_type::mass_type> ()>
wrap_get_mass(std::shared_ptr<particle_type const> self)
{
    return [=]() -> std::vector<typename particle_type::mass_type> {
        std::vector<typename particle_type::mass_type> output;
        output.reserve(self->nparticle());
        get_mass(*self, back_inserter(output));
        return std::move(output);
    };
}

template <typename particle_type>
static std::function<void (std::vector<typename particle_type::mass_type> const&)>
wrap_set_mass(std::shared_ptr<particle_type> self)
{
    return [=](std::vector<typename particle_type::mass_type> const& input) {
        if (input.size() != self->nparticle()) {
            throw std::invalid_argument("input array size not equal to number of particles");
        }
        set_mass(*self, input.begin());
    };
}

template <int dimension, typename float_type>
static int wrap_dimension(particle<dimension, float_type> const&)
{
    return dimension;
}

template <int dimension, typename float_type>
void particle<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    static std::string class_name = "particle_" + std::to_string(dimension);
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                class_<particle, std::shared_ptr<particle>>(class_name.c_str())
                    .def(constructor<size_type, unsigned int>())
                    .property("nparticle", &particle::nparticle)
                    .property("nspecies", &particle::nspecies)
                    .property("get_position", &wrap_get_position<particle>)
                    .property("set_position", &wrap_set_position<particle>)
                    .property("get_image", &wrap_get_image<particle>)
                    .property("set_image", &wrap_set_image<particle>)
                    .property("get_velocity", &wrap_get_velocity<particle>)
                    .property("set_velocity", &wrap_set_velocity<particle>)
                    .property("get_tag", &wrap_get_tag<particle>)
                    .property("set_tag", &wrap_set_tag<particle>)
                    .property("get_reverse_tag", &wrap_get_reverse_tag<particle>)
                    .property("set_reverse_tag", &wrap_set_reverse_tag<particle>)
                    .property("get_species", &wrap_get_species<particle>)
                    .property("set_species", &wrap_set_species<particle>)
                    .property("get_mass", &wrap_get_mass<particle>)
                    .property("set_mass", &wrap_set_mass<particle>)
                    .property("dimension", &wrap_dimension<dimension, float_type>)
                    .scope
                    [
                        class_<runtime>("runtime")
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

} // namespace gpu
} // namespace mdsim
} // namespace halmd
