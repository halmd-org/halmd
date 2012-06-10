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

#ifndef HALMD_MDSIM_GPU_PARTICLE_HPP
#define HALMD_MDSIM_GPU_PARTICLE_HPP

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/profiler.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <algorithm>
#include <vector>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension, typename float_type>
class particle
{
public:
    typedef typename type_traits<dimension, float_type>::vector_type vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type gpu_vector_type;

    typedef unsigned int size_type;
    typedef vector_type position_type;
    typedef vector_type image_type;
    typedef vector_type velocity_type;
    typedef unsigned int tag_type;
    typedef unsigned int reverse_tag_type;
    typedef unsigned int species_type;
    typedef float mass_type;

    typedef cuda::vector<float4> position_array_type;
    typedef cuda::vector<gpu_vector_type> image_array_type;
    typedef cuda::vector<float4> velocity_array_type;
    typedef cuda::vector<tag_type> tag_array_type;
    typedef cuda::vector<reverse_tag_type> reverse_tag_array_type;

    void rearrange(cuda::vector<unsigned int> const& g_index);

    /** grid and block dimensions for CUDA calls */
    cuda::config const dim;

    /**
     * Allocate particle arrays in GPU memory.
     *
     * @param nparticle number of particles
     * @param nspecies number of particle species
     *
     * All particle arrays, except the masses, are initialised to zero.
     * The particle masses are initialised to unit mass.
     */
    particle(size_type nparticle, unsigned int nspecies);

    /**
     * Returns number of particles.
     *
     * Currently the number of particles is fixed at construction of
     * particle. This may change in the future, to allow for chemical
     * reactions that do not conserve the number of particles, or to
     * transfer particles between domains of different processors.
     */
    size_type nparticle() const
    {
        return nparticle_;
    }

    /**
     * Returns number of species.
     */
    unsigned int nspecies() const
    {
        return nspecies_;
    }

    /**
     * Returns non-const reference to particle positions and species.
     */
    cache<position_array_type> const& position() const
    {
        return g_position_;
    }

    /**
     * Returns const reference to particle positions and species.
     */
    cache<position_array_type>& position()
    {
        return g_position_;
    }

    /**
     * Returns non-const reference to particle images.
     */
    cache<image_array_type> const& image() const
    {
        return g_image_;
    }

    /**
     * Returns const reference to particle images.
     */
    cache<image_array_type>& image()
    {
        return g_image_;
    }

    /**
     * Returns non-const reference to particle velocities and masses.
     */
    cache<velocity_array_type> const& velocity() const
    {
        return g_velocity_;
    }

    /**
     * Returns const reference to particle velocities and masses.
     */
    cache<velocity_array_type>& velocity()
    {
        return g_velocity_;
    }

    /**
     * Returns non-const reference to particle tags.
     */
    cache<tag_array_type> const& tag() const
    {
        return g_tag_;
    }

    /**
     * Returns const reference to particle tags.
     */
    cache<tag_array_type>& tag()
    {
        return g_tag_;
    }

    /**
     * Returns non-const reference to particle reverse tags.
     */
    cache<reverse_tag_array_type> const& reverse_tag() const
    {
        return g_reverse_tag_;
    }

    /**
     * Returns const reference to particle reverse tags.
     */
    cache<reverse_tag_array_type>& reverse_tag()
    {
        return g_reverse_tag_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** number of particles */
    size_type nparticle_;
    /** number of particle species */
    unsigned int nspecies_;
    /** positions, species */
    cache<position_array_type> g_position_;
    /** minimum image vectors */
    cache<image_array_type> g_image_;
    /** velocities, masses */
    cache<velocity_array_type> g_velocity_;
    /** particle tags */
    cache<tag_array_type> g_tag_;
    /** reverse particle tags */
    cache<reverse_tag_array_type> g_reverse_tag_;

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

/**
 * Copy particle positions to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_position(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::position_array_type position_array_type;
    cache_proxy<position_array_type const> g_position = particle.position();
    cuda::host::vector<typename position_array_type::value_type> h_position(g_position->size());
    cuda::copy(g_position->begin(), g_position->end(), h_position.begin());
    iterator_type output = first;
    for (typename position_array_type::value_type const& v : h_position) {
        typename particle_type::position_type position;
        typename particle_type::species_type species;
        tie(position, species) <<= v;
        *output++ = position;
    }
    return output;
}

/**
 * Copy particle positions from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_position(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::position_array_type position_array_type;
    cache_proxy<position_array_type> g_position = particle.position();
    cuda::host::vector<typename position_array_type::value_type> h_position(g_position->size());
    cuda::copy(g_position->begin(), g_position->end(), h_position.begin());
    iterator_type input = first;
    for (typename position_array_type::value_type& v : h_position) {
        typename particle_type::position_type position;
        typename particle_type::species_type species;
        tie(position, species) <<= v;
        position = *input++;
        v <<= tie(position, species);
    }
#ifdef USE_VERLET_DSFUN
    cuda::memset(g_position->begin(), g_position->begin() + g_position->capacity(), 0);
#endif
    cuda::copy(h_position.begin(), h_position.end(), g_position->begin());
    return input;
}

/**
 * Copy particle species to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_species(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::position_array_type position_array_type;
    cache_proxy<position_array_type const> g_position = particle.position();
    cuda::host::vector<typename position_array_type::value_type> h_position(g_position->size());
    cuda::copy(g_position->begin(), g_position->end(), h_position.begin());
    iterator_type output = first;
    for (typename position_array_type::value_type const& v : h_position) {
        typename particle_type::position_type position;
        typename particle_type::species_type species;
        tie(position, species) <<= v;
        *output++ = species;
    }
    return output;
}

/**
 * Copy particle species from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_species(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::position_array_type position_array_type;
    cache_proxy<position_array_type> g_position = particle.position();
    cuda::host::vector<typename position_array_type::value_type> h_position(g_position->size());
    cuda::copy(g_position->begin(), g_position->end(), h_position.begin());
    iterator_type input = first;
    for (typename position_array_type::value_type& v : h_position) {
        typename particle_type::position_type position;
        typename particle_type::species_type species;
        tie(position, species) <<= v;
        species = *input++;
        v <<= tie(position, species);
    }
    cuda::copy(h_position.begin(), h_position.end(), g_position->begin());
    return input;
}

/**
 * Copy particle images to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_image(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::image_array_type image_array_type;
    cache_proxy<image_array_type const> g_image = particle.image();
    cuda::host::vector<typename image_array_type::value_type> h_image(g_image->size());
    cuda::copy(g_image->begin(), g_image->end(), h_image.begin());
    return std::copy(h_image.begin(), h_image.end(), first);
}

/**
 * Copy particle images from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_image(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::image_array_type image_array_type;
    cache_proxy<image_array_type> g_image = particle.image();
    cuda::host::vector<typename image_array_type::value_type> h_image(g_image->size());
    iterator_type input = first;
    for (typename image_array_type::value_type& image : h_image) {
        image = *input++;
    }
    cuda::copy(h_image.begin(), h_image.end(), g_image->begin());
    return input;
}

/**
 * Copy particle velocities to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_velocity(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::velocity_array_type velocity_array_type;
    cache_proxy<velocity_array_type const> g_velocity = particle.velocity();
    cuda::host::vector<typename velocity_array_type::value_type> h_velocity(g_velocity->size());
    cuda::copy(g_velocity->begin(), g_velocity->end(), h_velocity.begin());
    iterator_type output = first;
    for (typename velocity_array_type::value_type const& v : h_velocity) {
        typename particle_type::velocity_type velocity;
        typename particle_type::mass_type mass;
        tie(velocity, mass) <<= v;
        *output++ = velocity;
    }
    return output;
}

/**
 * Copy particle velocities from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_velocity(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::velocity_array_type velocity_array_type;
    cache_proxy<velocity_array_type> g_velocity = particle.velocity();
    cuda::host::vector<typename velocity_array_type::value_type> h_velocity(g_velocity->size());
    cuda::copy(g_velocity->begin(), g_velocity->end(), h_velocity.begin());
    iterator_type input = first;
    for (typename velocity_array_type::value_type& v : h_velocity) {
        typename particle_type::velocity_type velocity;
        typename particle_type::mass_type mass;
        tie(velocity, mass) <<= v;
        velocity = *input++;
        v <<= tie(velocity, mass);
    }
#ifdef USE_VERLET_DSFUN
    cuda::memset(g_velocity->begin(), g_velocity->begin() + g_velocity->capacity(), 0);
#endif
    cuda::copy(h_velocity.begin(), h_velocity.end(), g_velocity->begin());
    return input;
}

/**
 * Copy particle masses to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_mass(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::velocity_array_type velocity_array_type;
    cache_proxy<velocity_array_type const> g_velocity = particle.velocity();
    cuda::host::vector<typename velocity_array_type::value_type> h_velocity(g_velocity->size());
    cuda::copy(g_velocity->begin(), g_velocity->end(), h_velocity.begin());
    iterator_type output = first;
    for (typename velocity_array_type::value_type const& v : h_velocity) {
        typename particle_type::velocity_type velocity;
        typename particle_type::mass_type mass;
        tie(velocity, mass) <<= v;
        *output++ = mass;
    }
    return output;
}

/**
 * Copy particle masses from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_mass(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::velocity_array_type velocity_array_type;
    cache_proxy<velocity_array_type> g_velocity = particle.velocity();
    cuda::host::vector<typename velocity_array_type::value_type> h_velocity(g_velocity->size());
    cuda::copy(g_velocity->begin(), g_velocity->end(), h_velocity.begin());
    iterator_type input = first;
    for (typename velocity_array_type::value_type& v : h_velocity) {
        typename particle_type::velocity_type velocity;
        typename particle_type::mass_type mass;
        tie(velocity, mass) <<= v;
        mass = *input++;
        v <<= tie(velocity, mass);
    }
    cuda::copy(h_velocity.begin(), h_velocity.end(), g_velocity->begin());
    return input;
}

/**
 * Copy particle tags to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_tag(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::tag_array_type tag_array_type;
    cache_proxy<tag_array_type const> g_tag = particle.tag();
    cuda::host::vector<typename tag_array_type::value_type> h_tag(g_tag->size());
    cuda::copy(g_tag->begin(), g_tag->end(), h_tag.begin());
    return std::copy(h_tag.begin(), h_tag.end(), first);
}

/**
 * Copy particle tags from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_tag(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::tag_array_type tag_array_type;
    cache_proxy<tag_array_type> g_tag = particle.tag();
    cuda::host::vector<typename tag_array_type::value_type> h_tag(g_tag->size());
    iterator_type input = first;
    for (typename tag_array_type::value_type& tag : h_tag) {
        tag = *input++;
    }
    cuda::copy(h_tag.begin(), h_tag.end(), g_tag->begin());
    return input;
}

/**
 * Copy particle reverse tags to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_reverse_tag(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::reverse_tag_array_type reverse_tag_array_type;
    cache_proxy<reverse_tag_array_type const> g_reverse_tag = particle.reverse_tag();
    cuda::host::vector<typename reverse_tag_array_type::value_type> h_reverse_tag(g_reverse_tag->size());
    cuda::copy(*g_reverse_tag, h_reverse_tag);
    return std::copy(h_reverse_tag.begin(), h_reverse_tag.end(), first);
}

/**
 * Copy particle reverse tags from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_reverse_tag(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::reverse_tag_array_type reverse_tag_array_type;
    cache_proxy<reverse_tag_array_type> g_reverse_tag = particle.reverse_tag();
    cuda::host::vector<typename reverse_tag_array_type::value_type> h_reverse_tag(g_reverse_tag->size());
    iterator_type input = first;
    for (typename reverse_tag_array_type::value_type& reverse_tag : h_reverse_tag) {
        reverse_tag = *input++;
    }
    cuda::copy(h_reverse_tag, *g_reverse_tag);
    return input;
}

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_HPP */
