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

#ifndef HALMD_MDSIM_HOST_PARTICLE_HPP
#define HALMD_MDSIM_HOST_PARTICLE_HPP

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/raw_array.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

#include <lua.hpp>

#include <algorithm>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class particle
{
public:
    typedef halmd::signal<void ()> signal_type;
    typedef signal_type::slot_function_type slot_function_type;

    typedef fixed_vector<float_type, dimension> vector_type;

    typedef unsigned int size_type;
    typedef vector_type position_type;
    typedef vector_type image_type;
    typedef vector_type velocity_type;
    typedef unsigned int tag_type;
    typedef unsigned int reverse_tag_type;
    typedef unsigned int species_type;
    typedef double mass_type;
    typedef fixed_vector<float_type, dimension> force_type;
    typedef double en_pot_type;
    typedef typename type_traits<dimension, float_type>::stress_tensor_type stress_pot_type;

    typedef raw_array<position_type> position_array_type;
    typedef raw_array<image_type> image_array_type;
    typedef raw_array<velocity_type> velocity_array_type;
    typedef raw_array<tag_type> tag_array_type;
    typedef raw_array<reverse_tag_type> reverse_tag_array_type;
    typedef raw_array<species_type> species_array_type;
    typedef raw_array<mass_type> mass_array_type;
    typedef raw_array<force_type> force_array_type;
    typedef raw_array<en_pot_type> en_pot_array_type;
    typedef raw_array<stress_pot_type> stress_pot_array_type;

    void rearrange(std::vector<unsigned int> const& index);

    /**
     * Allocate particle arrays in host memory.
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
     * Returns const reference to particle positions.
     */
    cache<position_array_type> const& position() const
    {
        return position_;
    }

    /**
     * Returns non-const reference to particle positions.
     */
    cache<position_array_type>& position()
    {
        return position_;
    }

    /**
     * Returns non-const reference to particle images.
     */
    cache<image_array_type> const& image() const
    {
        return image_;
    }

    /**
     * Returns const reference to particle images.
     */
    cache<image_array_type>& image()
    {
        return image_;
    }

    /**
     * Returns const reference to particle velocities.
     */
    cache<velocity_array_type> const& velocity() const
    {
        return velocity_;
    }

    /**
     * Returns non-const reference to particle velocities.
     */
    cache<velocity_array_type>& velocity()
    {
        return velocity_;
    }

    /**
     * Returns const reference to particle tags.
     */
    cache<tag_array_type> const& tag() const
    {
        return tag_;
    }

    /**
     * Returns non-const reference to particle tags.
     */
    cache<tag_array_type>& tag()
    {
        return tag_;
    }

    /**
     * Returns const reference to particle reverse tags.
     */
    cache<reverse_tag_array_type> const& reverse_tag() const
    {
        return reverse_tag_;
    }

    /**
     * Returns non-const reference to particle reverse tags.
     */
    cache<reverse_tag_array_type>& reverse_tag()
    {
        return reverse_tag_;
    }

    /**
     * Returns const reference to particle species.
     */
    cache<species_array_type> const& species() const
    {
        return species_;
    }

    /**
     * Returns non-const reference to particle species.
     */
    cache<species_array_type>& species()
    {
        return species_;
    }

    /**
     * Returns const reference to particle masses.
     */
    cache<mass_array_type> const& mass() const
    {
        return mass_;
    }

    /**
     * Returns non-const reference to particle masses.
     */
    cache<mass_array_type>& mass()
    {
        return mass_;
    }
        /**
     * Returns const reference to particle force.
     */
    cache<force_array_type> const& force()
    {
        update_force_();
        return force_;
    }

    /**
     * Returns non-const reference to particle force.
     */
    cache<force_array_type>& mutable_force()
    {
        return force_;
    }

    /**
     * Returns const reference to potential energy of particles.
     */
    cache<en_pot_array_type> const& potential_energy()
    {
        update_force_(true);
        return en_pot_;
    }

    /**
     * Returns non-const reference to potential energy of particles.
     */
    cache<en_pot_array_type>& mutable_potential_energy()
    {
        return en_pot_;
    }

    /**
     * Returns const reference to potential part of stress tensor.
     */
    cache<stress_pot_array_type> const& stress_pot()
    {
        update_force_(true);
        return stress_pot_;
    }

    /**
     * Returns non-const reference to potential part of stress tensor.
     */
    cache<stress_pot_array_type>& mutable_stress_pot()
    {
        return stress_pot_;
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
    unsigned int nparticle_;
    /** number of particle species */
    unsigned int nspecies_;
    /** positions, reduced to extended domain box */
    cache<position_array_type> position_;
    /** minimum image vectors */
    cache<image_array_type> image_;
    /** velocities */
    cache<velocity_array_type> velocity_;
    /** particle tags */
    cache<tag_array_type> tag_;
    /** reverse particle tags */
    cache<reverse_tag_array_type> reverse_tag_;
    /** particle species */
    cache<species_array_type> species_;
    /** particle masses */
    cache<mass_array_type> mass_;
    /** total force on particles */
    cache<force_array_type> force_;
    /** potential energy per particle */
    cache<en_pot_array_type> en_pot_;
    /** potential part of stress tensor for each particle */
    cache<stress_pot_array_type> stress_pot_;

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

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

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
    auto const& position = read_cache(particle.position());
    return std::copy(position.begin(), position.end(), first);
}

/**
 * Copy particle positions from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_position(particle_type& particle, iterator_type const& first)
{
    auto position = make_cache_mutable(particle.position());
    iterator_type input = first;
    for (auto& value : *position) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle images to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_image(particle_type const& particle, iterator_type const& first)
{
    auto const& image = read_cache(particle.image());
    return std::copy(image.begin(), image.end(), first);
}

/**
 * Copy particle images from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_image(particle_type& particle, iterator_type const& first)
{
    auto image = make_cache_mutable(particle.image());
    iterator_type input = first;
    for (auto& value : *image) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle velocities to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_velocity(particle_type const& particle, iterator_type const& first)
{
    auto const& velocity = read_cache(particle.velocity());
    return std::copy(velocity.begin(), velocity.end(), first);
}

/**
 * Copy particle velocities from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_velocity(particle_type& particle, iterator_type const& first)
{
    auto velocity = make_cache_mutable(particle.velocity());
    iterator_type input = first;
    for (auto& value : *velocity) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle tags to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_tag(particle_type const& particle, iterator_type const& first)
{
    auto const& tag = read_cache(particle.tag());
    return std::copy(tag.begin(), tag.end(), first);
}

/**
 * Copy particle tags from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_tag(particle_type& particle, iterator_type const& first)
{
    auto tag = make_cache_mutable(particle.tag());
    iterator_type input = first;
    for (auto& value : *tag) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle reverse tags to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_reverse_tag(particle_type const& particle, iterator_type const& first)
{
    auto const& reverse_tag = read_cache(particle.reverse_tag());
    return std::copy(reverse_tag.begin(), reverse_tag.end(), first);
}

/**
 * Copy particle reverse tags from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_reverse_tag(particle_type& particle, iterator_type const& first)
{
    auto reverse_tag = make_cache_mutable(particle.reverse_tag());
    iterator_type input = first;
    for (auto& value : *reverse_tag) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle species to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_species(particle_type const& particle, iterator_type const& first)
{
    auto const& species = read_cache(particle.species());
    return std::copy(species.begin(), species.end(), first);
}

/**
 * Copy particle species from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_species(particle_type& particle, iterator_type const& first)
{
    auto species = make_cache_mutable(particle.species());
    iterator_type input = first;
    for (auto& value : *species) {
        value = *input++;
    }
    return input;
}

/**
 * Copy particle masses to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_mass(particle_type const& particle, iterator_type const& first)
{
    auto const& mass = read_cache(particle.mass());
    return std::copy(mass.begin(), mass.end(), first);
}

/**
 * Copy particle masses from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_mass(particle_type& particle, iterator_type const& first)
{
    auto mass = make_cache_mutable(particle.mass());
    iterator_type input = first;
    for (auto& value : *mass) {
        value = *input++;
    }
    return input;
}

/**
 * Copy force per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_force(particle_type& particle, iterator_type const& first)
{
    return std::copy(particle.force()->begin(), particle.force()->end(), first);
}

/**
 * Copy potential energy per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_potential_energy(particle_type& particle, iterator_type const& first)
{
    return std::copy(particle.potential_energy()->begin(), particle.potential_energy()->end(), first);
}

/**
 * Copy potential part of stress tensor per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_stress_pot(particle_type& particle, iterator_type const& first)
{
    return std::copy(particle.stress_pot()->begin(), particle.stress_pot()->end(), first);
}


} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_HPP */
