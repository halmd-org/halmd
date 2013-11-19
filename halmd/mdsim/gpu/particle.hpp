/*
 * Copyright © 2010-2012 Felix Höfling
 * Copyright © 2013      Nicolas Höft
 * Copyright © 2008-2012 Peter Colberg
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
#include <halmd/utility/signal.hpp>
#include <halmd/mdsim/force_kernel.hpp>

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
    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

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
    typedef fixed_vector<float_type, dimension> force_type;
    typedef float_type en_pot_type;
    typedef typename type_traits<dimension, float_type>::stress_tensor_type stress_pot_type;

    typedef cuda::vector<float4> position_array_type;
    typedef cuda::vector<gpu_vector_type> image_array_type;
    typedef cuda::vector<float4> velocity_array_type;
    typedef cuda::vector<tag_type> tag_array_type;
    typedef cuda::vector<reverse_tag_type> reverse_tag_array_type;
    typedef cuda::vector<typename type_traits<dimension, float_type>::gpu::coalesced_vector_type> force_array_type;
    typedef cuda::vector<en_pot_type> en_pot_array_type;
    typedef cuda::vector<typename stress_pot_type::value_type> stress_pot_array_type;

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
     * Returns const reference to particle positions and species.
     */
    cache<position_array_type> const& position() const
    {
        return g_position_;
    }

    /**
     * Returns non-const reference to particle positions and species.
     */
    cache<position_array_type>& position()
    {
        return g_position_;
    }

    /**
     * Returns const reference to particle images.
     */
    cache<image_array_type> const& image() const
    {
        return g_image_;
    }

    /**
     * Returns non-const reference to particle images.
     */
    cache<image_array_type>& image()
    {
        return g_image_;
    }

    /**
     * Returns const reference to particle velocities and masses.
     */
    cache<velocity_array_type> const& velocity() const
    {
        return g_velocity_;
    }

    /**
     * Returns non-const reference to particle velocities and masses.
     */
    cache<velocity_array_type>& velocity()
    {
        return g_velocity_;
    }

    /**
     * Returns const reference to particle tags.
     */
    cache<tag_array_type> const& tag() const
    {
        return g_tag_;
    }

    /**
     * Returns non-const reference to particle tags.
     */
    cache<tag_array_type>& tag()
    {
        return g_tag_;
    }

    /**
     * Returns const reference to particle reverse tags.
     */
    cache<reverse_tag_array_type> const& reverse_tag() const
    {
        return g_reverse_tag_;
    }

    /**
     * Returns non-const reference to particle reverse tags.
     */
    cache<reverse_tag_array_type>& reverse_tag()
    {
        return g_reverse_tag_;
    }

    /**
     * Returns const reference to particle force.
     */
    cache<force_array_type> const& force()
    {
        update_force_();
        return g_force_;
    }

    /**
     * Returns non-const reference to particle force.
     */
    cache<force_array_type>& mutable_force()
    {
        return g_force_;
    }

    /**
     * Returns const reference to potential energy of particles.
     */
    cache<en_pot_array_type> const& potential_energy()
    {
        update_force_(true);
        return g_en_pot_;
    }

    /**
     * Returns non-const reference to potential energy of particles.
     */
    cache<en_pot_array_type>& mutable_potential_energy()
    {
        return g_en_pot_;
    }

    /**
     * Returns const reference to potential part of stress tensor.
     */
    cache<stress_pot_array_type> const& stress_pot()
    {
        update_force_(true);
        return g_stress_pot_;
    }
    /**
     * Returns non-const reference to potential part of stress tensor.
     */
    cache<stress_pot_array_type>& mutable_stress_pot()
    {
        return g_stress_pot_;
    }

    /**
     * Enable computation of auxiliary variables.
     *
     * The flag is reset after the next trigger of on_force_().
     */
    void aux_enable();

    /**
     * Returns true if computation of auxiliary variables is enabled.
     */
    bool aux_enabled() const
    {
        return aux_enabled_;
    }

    /**
     * Returns true if the force has to be reset to zero prior to reading.
     */
    bool force_zero()
    {
        return force_zero_;
    }

    /**
     * Disable a reset of the force to zero upon reading.
     *
     * Must be called after computation of the first force contribution.
     */
    void force_zero_disable()
    {
        force_zero_ = false;
    }

    /**
     * Indicate that a force update (ie. triggering on_force_()) is required
     */
    void mark_force_dirty()
    {
        force_dirty_ = true;

    }

    /**
     * Indicate that an auxiliary update (ie. triggering on_force_() with
     * aux_enabled) is required
     */
    void mark_aux_dirty()
    {
        aux_dirty_ = true;
    }

    connection on_prepend_force(slot_function_type const& slot)
    {
        return on_prepend_force_.connect(slot);
    }

    connection on_force(slot_function_type const& slot)
    {
        return on_force_.connect(slot);
    }

    connection on_append_force(slot_function_type const& slot)
    {
        return on_append_force_.connect(slot);
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
    /** total force on particles */
    cache<force_array_type> g_force_;
    /** potential energy per particle */
    cache<en_pot_array_type> g_en_pot_;
    /** potential part of stress tensor for each particle */
    cache<stress_pot_array_type> g_stress_pot_;

    /** flag that the force has to be reset to zero prior to reading */
    bool force_zero_;
    /** flag that the force cache is dirty (not up to date) */
    bool force_dirty_;
    /** flag that the caches of the auxiliary variables are dirty (not up to date) */
    bool aux_dirty_;
    /** flag that the computation of auxiliary variables is requested */
    bool aux_enabled_;

    /**
     * Update all forces and auxiliary variables if needed. The auxiliary
     * variables are guaranteed to be up-to-date upon return if with_aux was
     * set to true.
     *
     * Auxiliary variables are computed only if they are out of date
     * (aux_dirty_ == true) and if either with_aux or aux_enabled_ is true.
     *
     * Emit a warning if the force update would be necessary solely to compute
     * the auxiliary variables, which indicates a performance problem.
     */
    void update_force_(bool with_aux=false);

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type rearrange;
    };

    /** profiling runtime accumulators */
    runtime runtime_;

    signal_type on_prepend_force_;
    signal_type on_force_;
    signal_type on_append_force_;
};

/**
 * Copy particle positions to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_position(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::position_array_type position_array_type;
    position_array_type const& g_position = read_cache(particle.position());
    cuda::host::vector<typename position_array_type::value_type> h_position(g_position.size());
    cuda::copy(g_position.begin(), g_position.end(), h_position.begin());
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
    auto g_position = make_cache_mutable(particle.position());
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
    position_array_type const& g_position = read_cache(particle.position());
    cuda::host::vector<typename position_array_type::value_type> h_position(g_position.size());
    cuda::copy(g_position.begin(), g_position.end(), h_position.begin());
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
    auto g_position = make_cache_mutable(particle.position());
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
    image_array_type const& g_image = read_cache(particle.image());
    cuda::host::vector<typename image_array_type::value_type> h_image(g_image.size());
    cuda::copy(g_image.begin(), g_image.end(), h_image.begin());
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
    auto g_image = make_cache_mutable(particle.image());
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
    velocity_array_type const& g_velocity = read_cache(particle.velocity());
    cuda::host::vector<typename velocity_array_type::value_type> h_velocity(g_velocity.size());
    cuda::copy(g_velocity.begin(), g_velocity.end(), h_velocity.begin());
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
    auto g_velocity = make_cache_mutable(particle.velocity());
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
    velocity_array_type const& g_velocity = read_cache(particle.velocity());
    cuda::host::vector<typename velocity_array_type::value_type> h_velocity(g_velocity.size());
    cuda::copy(g_velocity.begin(), g_velocity.end(), h_velocity.begin());
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
    auto g_velocity = make_cache_mutable(particle.velocity());
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
    tag_array_type const& g_tag = read_cache(particle.tag());
    cuda::host::vector<typename tag_array_type::value_type> h_tag(g_tag.size());
    cuda::copy(g_tag.begin(), g_tag.end(), h_tag.begin());
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
    auto g_tag = make_cache_mutable(particle.tag());
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
    reverse_tag_array_type const& g_reverse_tag = read_cache(particle.reverse_tag());
    cuda::host::vector<typename reverse_tag_array_type::value_type> h_reverse_tag(g_reverse_tag.size());
    cuda::copy(g_reverse_tag.begin(), g_reverse_tag.end(), h_reverse_tag.begin());
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
    auto g_reverse_tag = make_cache_mutable(particle.reverse_tag());
    cuda::host::vector<typename reverse_tag_array_type::value_type> h_reverse_tag(g_reverse_tag->size());
    iterator_type input = first;
    for (typename reverse_tag_array_type::value_type& reverse_tag : h_reverse_tag) {
        reverse_tag = *input++;
    }
    cuda::copy(h_reverse_tag.begin(), h_reverse_tag.end(), g_reverse_tag->begin());
    return input;
}

/**
 * Copy net force per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_force(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::force_array_type force_array_type;
    force_array_type const& g_force = read_cache(particle.force());
    cuda::host::vector<typename force_array_type::value_type> h_force(g_force.size());
    cuda::copy(g_force.begin(), g_force.end(), h_force.begin());
    return std::copy(h_force.begin(), h_force.end(), first);
}

/**
 * Copy potential energy per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_potential_energy(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::en_pot_array_type en_pot_array_type;
    en_pot_array_type const& g_en_pot = read_cache(particle.potential_energy());
    cuda::host::vector<typename en_pot_array_type::value_type> h_en_pot(g_en_pot.size());
    cuda::copy(g_en_pot.begin(), g_en_pot.end(), h_en_pot.begin());
    return std::copy(h_en_pot.begin(), h_en_pot.end(), first);
}

/**
 * Copy potential part of stress tensor per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_stress_pot(particle_type& particle, iterator_type const& first)
{
    // copy data from GPU to host
    typedef typename particle_type::stress_pot_array_type stress_pot_array_type;
    stress_pot_array_type const& g_stress_pot = read_cache(particle.stress_pot());
    cuda::host::vector<typename stress_pot_array_type::value_type> h_stress_pot(g_stress_pot.size());
    h_stress_pot.reserve(g_stress_pot.capacity());
    cuda::copy(g_stress_pot.begin(), g_stress_pot.begin() + g_stress_pot.capacity(), h_stress_pot.begin());

    // convert from column-major to row-major layout
    typedef typename particle_type::stress_pot_type stress_pot_type;
    unsigned int stride = h_stress_pot.capacity() / stress_pot_type::static_size;
    iterator_type output = first;
    for (auto const& stress : h_stress_pot) {
        *output++ = read_stress_tensor<stress_pot_type>(&stress, stride);
    }
    return output;
}

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_HPP */
