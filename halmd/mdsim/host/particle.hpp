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

#ifndef HALMD_MDSIM_HOST_PARTICLE_HPP
#define HALMD_MDSIM_HOST_PARTICLE_HPP

#include <algorithm> // std::copy
#include <lua.hpp>

#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class particle
{
public:
    typedef fixed_vector<float_type, dimension> vector_type;

    typedef vector_type position_type;
    typedef vector_type image_type;
    typedef vector_type velocity_type;
    typedef unsigned int tag_type;
    typedef unsigned int reverse_tag_type;
    typedef unsigned int species_type;
    typedef double mass_type;
    typedef vector_type force_type;
    typedef double en_pot_type;
    typedef typename type_traits<dimension, float_type>::stress_tensor_type stress_pot_type;
    typedef double hypervirial_type;

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
    typedef raw_array<hypervirial_type> hypervirial_array_type;

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
    particle(std::size_t nparticle, unsigned int nspecies);

    /**
     * Returns number of particles.
     *
     * Currently the number of particles is fixed at construction of
     * particle. This may change in the future, to allow for chemical
     * reactions that do not conserve the number of particles, or to
     * transfer particles between domains of different processors.
     */
    std::size_t nparticle() const
    {
        return tag_.size();
    }

    /**
     * Returns number of particle placeholders.
     *
     * Currently the number of placeholders, i.e. the element count of the
     * particle arrays in memory, is equal to the number of particles.
     */
    std::size_t nplaceholder() const
    {
        return tag_.size();
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
    position_array_type const& position() const
    {
        return position_;
    }

    /**
     * Returns non-const reference to particle positions.
     */
    position_array_type& position()
    {
        return position_;
    }

    /**
     * Returns non-const reference to particle images.
     */
    image_array_type const& image() const
    {
        return image_;
    }

    /**
     * Returns const reference to particle images.
     */
    image_array_type& image()
    {
        return image_;
    }

    /**
     * Returns const reference to particle velocities.
     */
    velocity_array_type const& velocity() const
    {
        return velocity_;
    }

    /**
     * Returns non-const reference to particle velocities.
     */
    velocity_array_type& velocity()
    {
        return velocity_;
    }

    /**
     * Returns const reference to particle tags.
     */
    tag_array_type const& tag() const
    {
        return tag_;
    }

    /**
     * Returns non-const reference to particle tags.
     */
    tag_array_type& tag()
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
    species_array_type const& species() const
    {
        return species_;
    }

    /**
     * Returns non-const reference to particle species.
     */
    species_array_type& species()
    {
        return species_;
    }

    /**
     * Returns const reference to particle masses.
     */
    mass_array_type const& mass() const
    {
        return mass_;
    }

    /**
     * Returns non-const reference to particle masses.
     */
    mass_array_type& mass()
    {
        return mass_;
    }

    /**
     * Returns non-const reference to force per particle.
     */
    force_array_type const& force() const
    {
        return force_;
    }

    /**
     * Returns const reference to force per particle.
     */
    force_array_type& force()
    {
        return force_;
    }

    /**
     * Returns const reference to potential energy per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    en_pot_array_type const& en_pot() const
    {
        return en_pot_;
    }

    /**
     * Returns non-const reference to potential energy per particle.
     */
    en_pot_array_type& en_pot()
    {
        return en_pot_;
    }

    /**
     * Returns const reference to potential part of stress tensor per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    stress_pot_array_type const& stress_pot() const
    {
        return stress_pot_;
    }

    /**
     * Returns non-const reference to potential part of stress tensor per particle.
     */
    stress_pot_array_type& stress_pot()
    {
        return stress_pot_;
    }

    /**
     * Returns const reference to hypervirial per particle.
     *
     * This method checks that the computation of auxiliary variables was enabled.
     */
    hypervirial_array_type const& hypervirial() const
    {
        return hypervirial_;
    }

    /**
     * Returns non-const reference to hypervirial per particle.
     */
    hypervirial_array_type& hypervirial()
    {
        return hypervirial_;
    }

    /**
     * Enable computation of auxiliary variables.
     *
     * The flag is reset by the next call to prepare().
     */
    void aux_enable();

    /**
     * Returns true if computation of auxiliary variables is enabled.
     */
    bool aux_valid() const
    {
        return aux_valid_;
    }

    /**
     * Reset forces, and optionally auxiliary variables, to zero.
     */
    void prepare();

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** number of particle species */
    unsigned int nspecies_;
    /** positions, reduced to extended domain box */
    position_array_type position_;
    /** minimum image vectors */
    image_array_type image_;
    /** velocities */
    velocity_array_type velocity_;
    /** particle tags */
    tag_array_type tag_;
    /** reverse particle tags */
    cache<reverse_tag_array_type> reverse_tag_;
    /** particle species */
    species_array_type species_;
    /** particle masses */
    mass_array_type mass_;
    /** force per particle */
    force_array_type force_;
    /** potential energy per particle */
    en_pot_array_type en_pot_;
    /** potential part of stress tensor per particle */
    stress_pot_array_type stress_pot_;
    /** hypervirial per particle */
    hypervirial_array_type hypervirial_;

    /** flag for enabling the computation of auxiliary variables this step */
    bool aux_flag_;
    /** flag that indicates the auxiliary variables are computed this step */
    bool aux_valid_;

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
    position_array_type const& position = particle.position();
    return std::copy(position.begin(), position.end(), first);
}

/**
 * Copy particle positions from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_position(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::position_array_type position_array_type;
    position_array_type& position = particle.position();
    iterator_type input = first;
    for (typename particle_type::position_type& value : position) {
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
    typedef typename particle_type::image_array_type image_array_type;
    image_array_type const& image = particle.image();
    return std::copy(image.begin(), image.end(), first);
}

/**
 * Copy particle images from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_image(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::image_array_type image_array_type;
    image_array_type& image = particle.image();
    iterator_type input = first;
    for (typename particle_type::image_type& value : image) {
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
    typedef typename particle_type::velocity_array_type velocity_array_type;
    velocity_array_type const& velocity = particle.velocity();
    return std::copy(velocity.begin(), velocity.end(), first);
}

/**
 * Copy particle velocities from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_velocity(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::velocity_array_type velocity_array_type;
    velocity_array_type& velocity = particle.velocity();
    iterator_type input = first;
    for (typename particle_type::velocity_type& value : velocity) {
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
    typedef typename particle_type::tag_array_type tag_array_type;
    tag_array_type const& tag = particle.tag();
    return std::copy(tag.begin(), tag.end(), first);
}

/**
 * Copy particle tags from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_tag(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::tag_array_type tag_array_type;
    tag_array_type& tag = particle.tag();
    iterator_type input = first;
    for (typename particle_type::tag_type& value : tag) {
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
    typedef typename particle_type::reverse_tag_array_type reverse_tag_array_type;
    cache_proxy<reverse_tag_array_type const> reverse_tag = particle.reverse_tag();
    return std::copy(reverse_tag->begin(), reverse_tag->end(), first);
}

/**
 * Copy particle reverse tags from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_reverse_tag(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::reverse_tag_array_type reverse_tag_array_type;
    cache_proxy<reverse_tag_array_type> reverse_tag = particle.reverse_tag();
    iterator_type input = first;
    for (typename particle_type::reverse_tag_type& value : *reverse_tag) {
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
    typedef typename particle_type::species_array_type species_array_type;
    species_array_type const& species = particle.species();
    return std::copy(species.begin(), species.end(), first);
}

/**
 * Copy particle species from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_species(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::species_array_type species_array_type;
    species_array_type& species = particle.species();
    iterator_type input = first;
    for (typename particle_type::species_type& value : species) {
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
    typedef typename particle_type::mass_array_type mass_array_type;
    mass_array_type const& mass = particle.mass();
    return std::copy(mass.begin(), mass.end(), first);
}

/**
 * Copy particle masses from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_mass(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::mass_array_type mass_array_type;
    mass_array_type& mass = particle.mass();
    iterator_type input = first;
    for (typename particle_type::mass_type& value : mass) {
        value = *input++;
    }
    return input;
}

/**
 * Copy force per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_force(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::force_array_type force_array_type;
    force_array_type const& force = particle.force();
    return std::copy(force.begin(), force.end(), first);
}

/**
 * Copy force per particle from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_force(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::force_array_type force_array_type;
    force_array_type& force = particle.force();
    iterator_type input = first;
    for (typename particle_type::force_type& value : force) {
        value = *input++;
    }
    return input;
}

/**
 * Copy potential energy per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_en_pot(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::en_pot_array_type en_pot_array_type;
    en_pot_array_type const& en_pot = particle.en_pot();
    return std::copy(en_pot.begin(), en_pot.end(), first);
}

/**
 * Copy potential energy per particle from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_en_pot(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::en_pot_array_type en_pot_array_type;
    en_pot_array_type& en_pot = particle.en_pot();
    iterator_type input = first;
    for (typename particle_type::en_pot_type& value : en_pot) {
        value = *input++;
    }
    return input;
}

/**
 * Copy potential part of stress tensor per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_stress_pot(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::stress_pot_array_type stress_pot_array_type;
    stress_pot_array_type const& stress_pot = particle.stress_pot();
    return std::copy(stress_pot.begin(), stress_pot.end(), first);
}

/**
 * Copy potential part of stress tensor per particle from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_stress_pot(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::stress_pot_array_type stress_pot_array_type;
    stress_pot_array_type& stress_pot = particle.stress_pot();
    iterator_type input = first;
    for (typename particle_type::stress_pot_type& value : stress_pot) {
        value = *input++;
    }
    return input;
}

/**
 * Copy hypervirial per particle to given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
get_hypervirial(particle_type const& particle, iterator_type const& first)
{
    typedef typename particle_type::hypervirial_array_type hypervirial_array_type;
    hypervirial_array_type const& hypervirial = particle.hypervirial();
    return std::copy(hypervirial.begin(), hypervirial.end(), first);
}

/**
 * Copy hypervirial per particle from given array.
 */
template <typename particle_type, typename iterator_type>
inline iterator_type
set_hypervirial(particle_type& particle, iterator_type const& first)
{
    typedef typename particle_type::hypervirial_array_type hypervirial_array_type;
    hypervirial_array_type& hypervirial = particle.hypervirial();
    iterator_type input = first;
    for (typename particle_type::hypervirial_type& value : hypervirial) {
        value = *input++;
    }
    return input;
}

} // namespace mdsim
} // namespace host
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_HPP */
